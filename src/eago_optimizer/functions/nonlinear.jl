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


"""
$(FUNCTIONAME)

Stores a general quadratic function with a buffer.
"""
mutable struct BufferedNonlinear{V} <: AbstractEAGOConstraint
    "List of nodes in nonlinear expression"
    nd::Vector{JuMP.NodeData}
    "Adjacency Matrix for the expression"
    adj::SparseMatrixCSC{Bool, Int64}
    const_values::Vector{Float64}
    setstorage::Vector{V}
    numberstorage::Vector{Float64}
    isnumber::Vector{Bool}
    tp1storage::Vector{Float64}
    tp2storage::Vector{Float64}
    tp3storage::Vector{Float64}
    tp4storage::Vector{Float64}
    tpdict::Dict{Int64,Tuple{Int64,Int64,Int64,Int64}}
    grad_sparsity::Vector{Int64}
    dependent_variable_count::Int
    dependent_subexpressions::Vector{Int64}
    bnds::MOI.NLPBoundsPair
    linearity::JuMP._Derivatives.Linearity
    buffer::SAF
end

function BufferedNonlinear(func::JuMP._FunctionStorage, bnds::MOI.NLPBoundsPair, tag::T) where T <: RelaxTag

    nd = copy(func.nd)
    adj = copy(func.adj)

    const_values = copy(func.const_values)

    setstorage = fill(MC{N,T}(Interval(-Inf, Inf)), lenx)
    numberstorage = zeros(lenx)
    isnumber = fill(false, lenx)

    tpdict = Dict{Int64,Tuple{Int64,Int64,Int64,Int64}}()
    tp1_count = 0
    tp2_count = 0
    for i = 1:lenx
        node = @inbounds nd[i]
        op = node.index
        if double_tp(op)
            tp1_count += 1
            tpdict[i] = (tp1_count, tp1_count, tp2_count, tp2_count)
        elseif single_tp(op)
            tp1_count += 1
            tp2_count += 1
            tpdict[i] = (tp1_count, tp1_count, -1, -1)
        end
    end
    tp1storage = zeros(tp1_count)
    tp2storage = zeros(tp1_count)
    tp3storage = zeros(tp2_count)
    tp4storage = zeros(tp2_count)

    dependent_variable_count = length(func.grad_sparsity)
    grad_sparsity = copy(func.grad_sparsity)
    sort!(grad_sparsity)

    dependent_subexpressions = copy(func.dependent_subexpressions)

    saf_buffer = SAF(SAT[SAT(0.0, VI(-1)) for i=1:dependent_variable_count], 0.0)

    return BufferedNonlinear{MC{N,T}}(nd, adj, const_values, setstorage, numberstorage, isnumber,
                                      grad_sparsity, dependent_variable_count, dependent_subexpressions,
                                      bnds, saf_buffer)
end

function BufferedNonlinear{V}()
    return BufferedNonlinear{MC{N,T}}(JuMP.NodeData[], spzeros(Bool, 1), Float64[], setstorage,
                                      Float64[], Bool[], Int64[], 0, Int64[],
                                      MOI.NLPBoundsPair(-Inf, Inf), SAF(SAT[], 0.0))
end

"""
$(FUNCTIONAME)

Stores a general quadratic function with a buffer.
"""
mutable struct BufferedNonlinearSubexpression <: AbstractEAGOConstraint
    nd::Vector{JuMP.NodeData}
    adj::SparseMatrixCSC{Bool,Int64}
    const_values::Vector{Float64}
    setstorage::Vector{MC{N,T}}
    numberstorage::Vector{Float64}
    tp1storage::Vector{Float64}
    tp2storage::Vector{Float64}
    tp3storage::Vector{Float64}
    tp4storage::Vector{Float64}
    tpdict::Dict{Int64,Tuple{Int64,Int64,Int64,Int64}}
    linearity::JuMP._Derivatives.Linearity
end

function BufferedNonlinearSubexpression{V}(sub::JuMP._SubexpressionStorage)
end
