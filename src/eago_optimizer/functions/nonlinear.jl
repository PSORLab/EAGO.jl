# Copyright (c) 2018: Matthew Wilhelm & Matthew Stuber.
# This work is licensed under the Creative Commons Attribution-NonCommercial-
# ShareAlike 4.0 International License. To view a copy of this license, visit
# http://creativecommons.org/licenses/by-nc-sa/4.0/ or send a letter to Creative
# Commons, PO Box 1866, Mountain View, CA 94042, USA.
#############################################################################
# EAGO
# A development environment for robust and global optimization
# See https://github.com/PSORLab/EAGO.jl
#############################################################################
# TODO
#############################################################################

#=
NodeType convention is defined to parallel JuMP's nonlinear evaluator
@enum NodeType: CALL CALLUNIVAR VARIABLE VALUE SUBEXPRESSION PARAMETER

const OPERATORS = [:+, :-, :*, :^, :/, :ifelse, :max, :min]
const USER_OPERATOR_ID_START = length(operators) + 1
const OPERATOR_TO_ID = Dict{Symbol,Int}()
for i = 1:length(OPERATORS)
    OPERATOR_TO_ID[OPERATORS[i]] = i
end

const UNIVARIATE_OPERATORS = Symbol[:+, :-, :abs]
const UNIVARIATE_OPERATOR_TO_ID = Dict{Symbol,Int}(:+ => 1, :- => 2, :abs => 3)
const UNIVARIATE_OPERATOR_DERIV = Any[:(one(x)), :(-one(x)), :(ifelse(x >= 0, one(x), -one(x)))]
for (op, deriv) in Calculus.symbolic_derivatives_1arg()
    push!(UNIVARIATE_OPERATORS, op)
    push!(UNIVARIATE_OPERATOR_DERIV, deriv)
    UNIVARIATE_OPERATOR_TO_ID[op] = length(UNIVARIATE_OPERATORS)
end
const USER_UNIVAR_OPERATOR_ID_START = length(UNIVARIATE_OPERATORS) + 1

"""
$(TYPDEF)

Stores a general nonlinear function with a buffer.
"""
mutable struct BufferedNonlinearEq{N, T<:RelaxTag} <: AbstractEAGOConstraint
    "Main evaluator"
    evaluator::Evaluator
    "List of nodes in nonlinear expression"
    node_list::NodeData
    const_values::Vector{Float64}
    set_storage::Vector{MC{N,T}}
    grad_sparsity::Vector{Int64}
    dependent_subexpressions::Vector{Int64}
end

function interval_bound(f::BufferedNonlinearEq{N,T}, n::NodeBB) where {N, T<:RelaxTag}
end

function relax!(m::Optimizer, f::BufferedNonlinearEq{N,T}, indx::Int, check_safe::Bool) where {N, T<:RelaxTag}
end

"""
$(FUNCTIONAME)

Stores a general quadratic function with a buffer.
"""
mutable struct BufferedNonlinear <: AbstractEAGOConstraint
    func
end
=#

####
#### Nonlinear Storage
####

struct EmptyNLPEvaluator <: MOI.AbstractNLPEvaluator
    _current_node::NodeBB
    has_nlobj::Bool
end
EmptyNLPEvaluator() = EmptyNLPEvaluator(NodeBB(),false)
set_current_node!(x::EmptyNLPEvaluator, n::NodeBB) = ()

MOI.features_available(::EmptyNLPEvaluator) = [:Grad, :Jac, :Hess]
MOI.initialize(::EmptyNLPEvaluator, features) = nothing
MOI.eval_objective(::EmptyNLPEvaluator, x) = NaN
function MOI.eval_constraint(::EmptyNLPEvaluator, g, x)
    @assert length(g) == 0
    return
end
MOI.eval_objective_gradient(::EmptyNLPEvaluator, g, x) = nothing
MOI.jacobian_structure(::EmptyNLPEvaluator) = Tuple{Int64,Int64}[]
MOI.hessian_lagrangian_structure(::EmptyNLPEvaluator) = Tuple{Int64,Int64}[]
function MOI.eval_constraint_jacobian(::EmptyNLPEvaluator, J, x)
    @assert length(J) == 0
    return
end
function MOI.eval_hessian_lagrangian(::EmptyNLPEvaluator, H, x, σ, μ)
    @assert length(H) == 0
    return
end

empty_nlp_data() = MOI.NLPBlockData([], EmptyNLPEvaluator(), false)
