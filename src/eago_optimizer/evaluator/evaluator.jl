mutable struct FunctionSetStorage{N,T<:RelaxTag}
    nd::Vector{JuMP.NodeData}
    adj::SparseMatrixCSC{Bool,Int64}
    const_values::Vector{Float64}
    setstorage::Vector{MC{N,T}}
    numberstorage::Vector{Float64}
    numvalued::Vector{Bool}
    tp1storage::Vector{Float64}
    tp2storage::Vector{Float64}
    tp3storage::Vector{Float64}
    tp4storage::Vector{Float64}
    tpdict::Dict{Int64,Tuple{Int64,Int64,Int64,Int64}}
    grad_sparsity::Vector{Int64}
    hess_I::Vector{Int64}
    hess_J::Vector{Int64}
    dependent_subexpressions::Vector{Int64}
end

FunctionSetStorage(N,T) = FunctionSetStorage{N,T}(JuMP.NodeData[],spzeros(Bool,1),
                                           Float64[],MC{N,T}[],Float64[], Bool[],
                                           Float64[], Float64[], Float64[], Float64[],
                                           Dict{Int64,Tuple{Int64,Int64,Int64,Int64}}(), Int[],Int[],Int[],Int[])

mutable struct SubexpressionSetStorage{N,T<:RelaxTag}
    nd::Vector{JuMP.NodeData}
    adj::SparseMatrixCSC{Bool,Int64}
    const_values::Vector{Float64}
    setstorage::Vector{MC{N,T}}
    numberstorage::Vector{Float64}
    numvalued::Vector{Bool}
    tp1storage::Vector{Float64}
    tp2storage::Vector{Float64}
    tp3storage::Vector{Float64}
    tp4storage::Vector{Float64}
    tpdict::Dict{Int64,Tuple{Int64,Int64,Int64,Int64}}
    linearity::JuMP._Derivatives.Linearity
end

function SubexpressionSetStorage(N::Int64, s::T, nd::Vector{JuMP.NodeData},
                                 const_values, num_variables,
                                 subexpression_linearity,
                                 moi_index_to_consecutive_index) where T<:RelaxTag

    nd = JuMP.replace_moi_variables(nd, moi_index_to_consecutive_index)
    len_nd = length(nd)
    adj = adjmat(nd)
    setstorage = zeros(MC{N}, len_nd)
    numberstorage = zeros(len_nd)
    numvalued = zeros(Bool, len_nd)
    linearity = JuMP.classify_linearity(nd, adj, subexpression_linearity)

    tpdict = Dict{Int64,Tuple{Int64,Int64,Int64,Int64}}()
    tp1_count = 0
    tp2_count = 0
    tp3_count = 0
    tp4_count = 0
    for i in 1:len_nd
        op = nd[i].index
        if double_tp(op)
            tp1_count += 1
            tp2_count += 1
            tp3_count += 1
            tp4_count += 1
            tpdict[i] = (tp1_count, tp2_count, tp3_count, tp4_count)
        elseif single_tp(op)
            tp1_count += 1
            tp2_count += 1
            tpdict[i] = (tp1_count, tp2_count, -1, -1)
        end
    end
    tp1storage = zeros(tp1_count)
    tp2storage = zeros(tp2_count)
    tp3storage = zeros(tp3_count)
    tp4storage = zeros(tp4_count)

    return SubexpressionSetStorage{N,T}(nd, adj, const_values, setstorage, numberstorage,
                                      numvalued, tp1storage, tp2storage, tp3storage,
                                      tp4storage, tpdict, linearity[1])
end

"""
    Evaluator

MOI.AbstractNLPEvaluator for calculating relaxations of nonlinear terms.
"""
mutable struct Evaluator{N, T<:RelaxTag} <: MOI.AbstractNLPEvaluator
    user_operators::JuMP._Derivatives.UserOperatorRegistry
    has_user_mv_operator::Bool
    parameter_values::Vector{Float64}
    variable_number::Int64
    index_to_variable::Vector{Tuple{Int64,Int64,Int64}}
    current_node::NodeBB
    disable_1storder::Bool
    constraint_number::Int64
    subexpression_number::Int64
    has_nlobj::Bool
    has_reverse::Bool
    subgrad_tighten::Bool
    subgrad_tighten_reverse::Bool
    first_eval_flag::Bool
    cp_reptitions::Int64
    cp_tolerance::Float64
    objective::FunctionSetStorage{N,T}
    objective_ubd::Float64
    constraints::Vector{FunctionSetStorage{N,T}}
    constraints_lbd::Vector{Float64}
    constraints_ubd::Vector{Float64}
    subexpressions::Vector{SubexpressionSetStorage{N,T}}
    subexpression_isnum::Vector{Bool}
    subexpression_order::Vector{Int64}
    subexpression_values_flt::Vector{Float64}
    subexpression_values_set::Vector{MC{N,T}}
    subexpression_linearity::Vector{JuMP._Derivatives.Linearity}
    last_x::Vector{Float64}
    last_node::NodeBB
    last_obj::MC{N,T}
    jac_storage::Vector{MC{N,T}}
    user_output_buffer::Vector{MC{N,T}}
    seeds::Vector{SVector{N,Float64}}
    function Evaluator{N,T}() where {N,T<:RelaxTag}
        d = new()
        d.user_operators = JuMP._Derivatives.UserOperatorRegistry()
        d.first_eval_flag = false
        d.objective_ubd = Inf
        d.constraints = FunctionSetStorage{N,T}[]
        d.constraints_lbd = Float64[]
        d.constraints_ubd = Float64[]
        d.objective = FunctionSetStorage(N,T)
        d.index_to_variable = Tuple{Int64,Int64,Int64}[]
        d.seeds = SVector{N,Float64}[]
        for i in 1:N
            push!(d.seeds, seed_gradient(i, N))
        end
        return d
    end
end

"""
    set_current_node!

Sets the current node in the Evaluator structure.
"""
function set_current_node!(x::Evaluator, n::NodeBB)
    x.current_node = n
end

include("univariate.jl")
include("passes.jl")
include("get_info.jl")
include("load.jl")
include("reform.jl")

num_state_variables(x::Evaluator) = 0
num_decision_variables(x::Evaluator) = x.variable_number
