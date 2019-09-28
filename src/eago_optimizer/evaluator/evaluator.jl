mutable struct FunctionSetStorage{N}
    nd::Vector{JuMP.NodeData}
    adj::SparseMatrixCSC{Bool,Int64}
    const_values::Vector{Float64}
    setstorage::Vector{MC{N}}
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

FunctionSetStorage(N) = FunctionSetStorage{N}(JuMP.NodeData[],spzeros(Bool,1),
                                           Float64[],MC{N}[],Float64[], Bool[],
                                           Float64[], Float64[], Float64[], Float64[],
                                           Dict{Int,Tuple{Int,Int,Int,Int}}(), Int[],Int[],Int[],Int[])
eltype(x::FunctionSetStorage{N}) where N = N

mutable struct SubexpressionSetStorage{N}
    nd::Vector{JuMP.NodeData}
    adj::SparseMatrixCSC{Bool,Int64}
    const_values::Vector{Float64}
    setstorage::Vector{MC{N}}
    numberstorage::Vector{Float64}
    numvalued::Vector{Bool}
    tp1storage::Vector{Float64}
    tp2storage::Vector{Float64}
    tp3storage::Vector{Float64}
    tp4storage::Vector{Float64}
    tpdict::Dict{Int64,Tuple{Int64,Int64,Int64,Int64}}
    linearity::JuMP._Derivatives.Linearity
end
eltype(x::SubexpressionSetStorage{N}) where N = N

function SubexpressionSetStorage(N::Int, nd::Vector{JuMP.NodeData},
                                 const_values, num_variables,
                                 subexpression_linearity,
                                 moi_index_to_consecutive_index)

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

    return SubexpressionSetStorage{N}(nd, adj, const_values, setstorage, numberstorage,
                                      numvalued, tp1storage, tp2storage, tp3storage,
                                      tp4storage, tpdict, linearity[1])
end

"""
    Evaluator

MOI.AbstractNLPEvaluator for calculating relaxations of nonlinear terms.
"""
mutable struct Evaluator{N} <: MOI.AbstractNLPEvaluator
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
    objective::FunctionSetStorage{N}
    objective_ubd::Float64
    constraints::Vector{FunctionSetStorage{N}}
    constraints_lbd::Vector{Float64}
    constraints_ubd::Vector{Float64}
    subexpressions::Vector{SubexpressionSetStorage{N}}
    subexpression_isnum::Vector{Bool}
    subexpression_order::Vector{Int64}
    subexpression_values_flt::Vector{Float64}
    subexpression_values_set::Vector{MC{N}}
    subexpression_linearity::Vector{JuMP._Derivatives.Linearity}
    last_x::Vector{Float64}
    last_node::NodeBB
    last_obj::MC{N}
    jac_storage::Vector{MC{N}}
    user_output_buffer::Vector{MC{N}}
    function Evaluator{N}() where N
        d = new()
        d.user_operators = JuMP._Derivatives.UserOperatorRegistry()
        d.first_eval_flag = false
        d.objective_ubd = Inf
        d.constraints = FunctionSetStorage{N}[]
        d.constraints_lbd = Float64[]
        d.constraints_ubd = Float64[]
        d.objective = FunctionSetStorage(N)
        d.index_to_variable = Tuple{Int64,Int64,Int64}[]
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
