mutable struct FunctionSetStorage{T}
    nd::Vector{JuMP.NodeData}
    adj::SparseMatrixCSC{Bool,Int}
    const_values::Vector{Float64}
    setstorage::Vector{T}
    numberstorage::Vector{Float64}
    numvalued::Vector{Bool}
    grad_sparsity::Vector{Int}
    hess_I::Vector{Int}
    hess_J::Vector{Int}
    dependent_subexpressions::Vector{Int}
end

FunctionSetStorage(T) = FunctionSetStorage{T}(JuMP.NodeData[],spzeros(Bool,1),
                                           Float64[],T[],Float64[],
                                           Bool[],Int[],Int[],Int[],Int[])

mutable struct SubexpressionSetStorage{T}
    nd::Vector{JuMP.NodeData}
    adj::SparseMatrixCSC{Bool,Int}
    const_values::Vector{Float64}
    setstorage::Vector{T}
    numberstorage::Vector{Float64}
    numvalued::Vector{Bool}
    linearity::JuMP._Derivatives.Linearity
end

function SubexpressionSetStorage(T,nd::Vector{JuMP.NodeData}, const_values, num_variables, subexpression_linearity, moi_index_to_consecutive_index)

    nd = JuMP.replace_moi_variables(nd, moi_index_to_consecutive_index)
    len_nd = length(nd)
    adj = adjmat(nd)
    setstorage = zeros(T,len_nd)
    numberstorage = zeros(len_nd)
    numvalued = zeros(Bool,len_nd)
    linearity = JuMP.classify_linearity(nd, adj, subexpression_linearity)

    return SubexpressionSetStorage{T}(nd, adj, const_values, setstorage, numberstorage,
                                      numvalued, linearity[1])
end

"""
    Evaluator

MOI.AbstractNLPEvaluator for calculating relaxations of nonlinear terms.
"""
mutable struct Evaluator{T<:Real} <: MOI.AbstractNLPEvaluator
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
    constraints::Vector{FunctionSetStorage{T}}
    constraints_lbd::Vector{Float64}
    constraints_ubd::Vector{Float64}
    subexpressions::Vector{SubexpressionSetStorage{T}}
    subexpression_isnum::Vector{Bool}
    subexpression_order::Vector{Int}
    subexpression_values_flt::Vector{Float64}
    subexpression_values_set::Vector{T}
    subexpression_linearity::Vector{JuMP._Derivatives.Linearity}
    subexpressions_as_julia_expressions::Vector{Any}
    last_x::Vector{Float64}
    last_node::NodeBB
    last_obj::T
    jac_storage::Vector{T}
    user_output_buffer::Vector{T}
    eval_objective_timer::Float64
    eval_constraint_timer::Float64
    eval_objective_gradient_timer::Float64
    eval_constraint_jacobian_timer::Float64
    eval_hessian_lagrangian_timer::Float64
    function Evaluator{T}(m) where T<:Real
        d = new()
        d.m = m
        d.first_eval_flag = false
        d.objective_ubd = Inf
        d.constraints = FunctionSetStorage[]
        d.constraints_lbd = Float64[]
        d.constraints_ubd = Float64[]
        d.objective = FunctionSetStorage(T)
        d.index_to_variable = Tuple{Int,Int,Int}[]
        d.eval_objective_timer = 0.0
        d.eval_constraint_timer = 0.0
        d.eval_objective_gradient_timer = 0.0
        d.eval_constraint_jacobian_timer = 0.0
        d.eval_hessian_lagrangian_timer = 0.0
        return d
    end
end

"""
    set_current_node!

Sets the current node in the Evaluator structure.
"""
function set_current_node!(x::Evaluator,n::NodeBB)
    x.current_node = n
end

eltype(x::Evaluator{T}) where T  = T


include("univariate.jl")
include("passes.jl")
include("get_info.jl")
include("load.jl")
include("reform.jl")

num_state_variables(x::Evaluator) = 0
num_decision_variables(x::Evaluator) = x.variable_number
