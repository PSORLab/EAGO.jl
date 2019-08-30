mutable struct FunctionSetStorage{N}
    nd::Vector{JuMP.NodeData}
    adj::SparseMatrixCSC{Bool,Int}
    const_values::Vector{Float64}
    setstorage::Vector{MC{N}}
    numberstorage::Vector{Float64}
    numvalued::Vector{Bool}
    tp1storage::Vector{Float64}
    tp2storage::Vector{Float64}
    tp3storage::Vector{Float64}
    tp4storage::Vector{Float64}
    tpdict::Dict{Int,Tuple{Int,Int,Int,Int}}
    grad_sparsity::Vector{Int}
    hess_I::Vector{Int}
    hess_J::Vector{Int}
    dependent_subexpressions::Vector{Int}
end

FunctionSetStorage(N) = FunctionSetStorage{N}(JuMP.NodeData[],spzeros(Bool,1),
                                           Float64[],T[],Float64[], Bool[],
                                           Float64[], Float64[], Float64[], Float64[],
                                           Int[],Int[],Int[],Int[])
eltype(x::FunctionSetStorage{N}) where N = N

mutable struct SubexpressionSetStorage{N}
    nd::Vector{JuMP.NodeData}
    adj::SparseMatrixCSC{Bool,Int}
    const_values::Vector{Float64}
    setstorage::Vector{MC{N}}
    numberstorage::Vector{Float64}
    numvalued::Vector{Bool}
    tp1storage::Vector{Float64}
    tp2storage::Vector{Float64}
    tp3storage::Vector{Float64}
    tp4storage::Vector{Float64}
    tpdict::Dict{Int,Tuple{Int,Int,Int,Int}}
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

    tpdict = Dict{Int,Tuple{Int,Int,Int,Int}}()
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
    m::Model
    has_user_mv_operator::Bool
    parameter_values::Vector{Float64}
    variable_number::Int
    index_to_variable::Vector{Tuple{Int,Int,Int}}
    current_node::NodeBB
    disable_1storder::Bool
    disable_2ndorder::Bool
    constraint_number::Int
    subexpression_number::Int
    has_nlobj::Bool
    has_reverse::Bool
    subgrad_tighten::Bool
    subgrad_tighten_reverse::Bool
    first_eval_flag::Bool
    cp_reptitions::Int
    cp_tolerance::Float64
    objective::FunctionSetStorage
    objective_ubd::Float64
    constraints::Vector{FunctionSetStorage{N}}
    constraints_lbd::Vector{Float64}
    constraints_ubd::Vector{Float64}
    subexpressions::Vector{SubexpressionSetStorage{N}}
    subexpression_isnum::Vector{Bool}
    subexpression_order::Vector{Int}
    subexpression_values_flt::Vector{Float64}
    subexpression_values_set::Vector{MC{N}}
    subexpression_linearity::Vector{JuMP._Derivatives.Linearity}
    subexpressions_as_julia_expressions::Vector{Any}
    last_x::Vector{Float64}
    last_node::NodeBB
    last_obj::MC{N}
    jac_storage::Vector{MC{N}}
    user_output_buffer::Vector{MC{N}}
    eval_objective_timer::Float64
    eval_constraint_timer::Float64
    eval_objective_gradient_timer::Float64
    eval_constraint_jacobian_timer::Float64
    eval_hessian_lagrangian_timer::Float64
    function Evaluator{N}(m) where N
        d = new()
        d.m = m
        d.first_eval_flag = false
        d.objective_ubd = Inf
        d.constraints = FunctionSetStorage[]
        d.constraints_lbd = Float64[]
        d.constraints_ubd = Float64[]
        d.objective = FunctionSetStorage(N)
        d.index_to_variable = Tuple{Int,Int,Int}[]
        d.eval_objective_timer = 0.0
        d.eval_constraint_timer = 0.0
        d.eval_objective_gradient_timer = 0.0
        d.eval_constraint_jacobian_timer = 0.0
        d.eval_hessian_lagrangian_timer = 0.0
        return d
    end
end
eltype(x::Evaluator{N}) where N  = N

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
