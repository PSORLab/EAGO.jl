# Copyright (c) 2018: Matthew Wilhelm & Matthew Stuber.
# This code is licensed under MIT license (see LICENSE.md for full details)
#############################################################################
# EAGO
# A development environment for robust and global optimization
# See https://github.com/PSORLab/EAGO.jl
#############################################################################
# Defines the NonlinearExpression, BufferedNonlinearFunction used in
# constructing relaxations of nonlinear functions along with a number of
# helper functions including an Evaluator structure and: set_node_flag!
# set_node!, set_reference_point!, retrieve_node, prior_eval
# copy_subexpression_value!, eliminate_fixed_variables!
#############################################################################

include(joinpath(@__DIR__, "empty_evaluator.jl"))
include(joinpath(@__DIR__, "register_special.jl"))
include(joinpath(@__DIR__, "graph", "abstract_graph.jl"))
include(joinpath(@__DIR__, "composite_relax", "composite_relax.jl"))

function linearity(d::JuMP._Derivatives.Linearity)
    if d == JuMP._Derivatives.LINEAR
        return Graph.LIN_LINEAR
    elseif d == JuMP._Derivatives.PIECEWISE_LINEAR
        return Graph.LIN_PIECEWISE_LINEAR
    elseif d == JuMP._Derivatives.NONLINEAR
        return Graph.LIN_NONLINEAR
    end
    return Graph.LIN_CONSTANT      # if d == JuMP._Derivatives.CONSTANT
end
function linearity(nd::Vector{JuMP._Derivatives.NodeData},
                   adj::SparseMatrixCSC{Int,Int},
                   d::Vector{JuMP._Derivatives.Linearity})
    x = JuMP._Derivatives.classify_linearity(nd, adj, d)
    linearity(x)
end

"""
$(TYPEDEF)

Stores a general quadratic function with a buffer.
"""
mutable struct NonlinearExpression{N,T} <: AbstractEAGOConstraint
    g::DirectedAcyclicGraph{Float64}
    relax_cache::RelaxCache{N,T}
end
function NonlinearExpression()
    g = DirectedAcyclicGraph{Float64}()
    c = RelaxCache{1,NS}()
    return NonlinearExpression{MC{1,NS}}(g, c)
end

"""
$(TYPEDEF)

Stores a general nonlinear function with a buffer represented by the sum of a tape
and a scalar affine function.
"""
mutable struct BufferedNonlinearFunction{N,T} <: AbstractEAGOConstraint
    expr::NonlinearExpression{N,T}
    saf::SAF
    lower_bound::Float64
    upper_bound::Float64
end
function BufferedNonlinearFunction()
    ex = NonlinearExpression()
    saf = SAF(SAT[], 0.0)
    return BufferedNonlinearFunction{1,NS}(ex, saf, 0.0, 0.0)
end

function NonlinearExpression!(sub::Union{JuMP._SubexpressionStorage,JuMP._FunctionStorage},
                              sub_sparsity::Dict{Int,Vector{Int}}, subexpr_indx::Int,
                              subexpr_linearity::Vector{JuMP._Derivatives.Linearity},
                              tag::T; is_sub::Bool = false) where T
    g = DirectedAcyclicGraph{Float64}(sub, sub_sparsity, subexpr_linearity)
    grad_sparsity = _sparsity(g)
    n = length(grad_sparsity)
    if is_sub
        sub_sparsity[subexpr_indx] = copy(grad_sparsity) # updates subexpression sparsity dictionary
    end
    c = RelaxCache{n,T}()
    return NonlinearExpression{MC{n,T}}(g, c)
end
function BufferedNonlinearFunction(f::JuMP._FunctionStorage, b::MOI.NLPBoundsPair,
                                   sub_sparsity::Dict{Int,Vector{Int}},
                                   subexpr_lin::Vector{JuMP._Derivatives.Linearity},
                                   tag::T) where T <: RelaxTag

    ex = NonlinearExpression!(f, sub_sparsity, subexpr_indx, subexpr_lin, tag)
    n = length(_sparsity(ex.g))
    saf = SAF(SAT[SAT(0.0, VI(i)) for i = 1:n], 0.0)
    return BufferedNonlinearFunction{MC{n,T}}(ex, saf, b.lower, b.upper)
end

function set_intersect_value!(expr::NonlinearExpression{V}, value) where V
    if !expr.isnumber[1]
        expr.value = expr.setstorage[1] ∩ value
        expr.setstorage[1] = expr.value
    end
    return
end

"""
    Evaluator

MOI.AbstractNLPEvaluator for calculating relaxations of nonlinear terms.

Checks that the resulting value should be a number...

$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct Evaluator <: MOI.AbstractNLPEvaluator
    user_operators::JuMP._Derivatives.UserOperatorRegistry = JuMP._Derivatives.UserOperatorRegistry()
    has_user_mv_operator::Bool = false
    num_mv_buffer::Vector{Float64} = Float64[]
    parameter_values::Vector{Float64} = Float64[]
    node::NodeBB = NodeBB()
    variable_values::VariableValues{Float64}
    subgrad_tighten::Bool = false
    reverse_subgrad_tighten::Bool = false
    ctx::GuardCtx = GuardCtx()
    subexpressions::Vector{NonlinearExpression} = NonlinearExpression[]
    subexpressions_eval::Vector{Bool}           = Bool[]
    is_post::Bool = false
    is_intersect::Bool = false
    is_first_eval::Bool = false
    interval_intersect::Bool = false
    subgrad_tol::Float64 = 1E-10
end

"""
$(FUNCTIONNAME)

Sets the current node in the Evaluator structure.
"""
function set_node!(d::Evaluator, n::NodeBB)
    d.node = NodeBB(n)
    @inbounds for i = 1:length(n)
        vi = d.node_to_variable_map[i]
        d.lower_variable_bounds[vi] = n.lower_variable_bounds[i]
        d.upper_variable_bounds[vi] = n.upper_variable_bounds[i]
    end
    fill!(d.subexpressions_eval, false)
    d.is_first_eval = true
    return nothing
end

function retrieve_node(d::Evaluator)
    n = d.current_node
    nv_map = d.node_to_variable_map
    return NodeBB(copy(d.lower_variable_bounds[nv_map]),
                  copy(d.upper_variable_bounds[nv_map]),
                  n.lower_bound, n.upper_bound, n.depth, n.id)
end
function retrieve_x!(out::Vector{Float64}, d::Evaluator)
    @inbounds for i = 1:length(d.node_to_variable_map)
        out[i] = d.x[d.node_to_variable_map[i]]
    end
    return nothing
end
prior_eval(d::Evaluator, i::Int) = d.subexpressions_eval[i]

#=
Assumes the sparsities are sorted...
=#
function copy_subexpression_value!(k::Int, op::Int, subexpression::NonlinearExpression{MC{N1,T}},
                                   numvalued::Vector{Bool}, numberstorage::Vector{S}, setstorage::Vector{MC{N2,T}},
                                   cv_buffer::Vector{S}, cc_buffer::Vector{S},
                                   func_sparsity::Vector{Int}) where {N1, N2, S, T <: RelaxTag}

    # fill cv_grad/cc_grad buffers
    sub_sparsity = subexpression.grad_sparsity
    sset = subexpression.setstorage[1]
    fill!(cv_buffer, zero(S))
    fill!(cc_buffer, zero(S))

    sub_sparsity_count = 1
    subs_index = @inbounds sub_sparsity[1]
    for i = 1:length(func_sparsity)
        func_index = @inbounds func_sparsity[i]
        if func_index === subs_index
            @inbounds cv_buffer[i] = sset.cv_grad[sub_sparsity_count]
            @inbounds cc_buffer[i] = sset.cc_grad[sub_sparsity_count]
            sub_sparsity_count += 1
            subs_index = @inbounds sub_sparsity[sub_sparsity_count]
        end
    end

    cv_grad = SVector(cv_buffer)
    cc_grad = SVector(cc_buffer)

    setstorage[k] = MC{N1,T}(sset.cv, sset.cc, sset.Intv, cv_grad, cc_grad, sset.cnst)

    return nothing
end

function eliminate_fixed_variables!(f::NonlinearExpression{V}, v::Vector{VariableInfo}) where V
    num_constants = length(f.const_values)
    indx_to_const_loc = Dict{Int,Int}()
    for i = 1:length(expr.nd)
        nd = @inbounds expr.nd[i]
        if nd.nodetype === JuMP._Derivatives.VARIABLE
            indx = nd.index
            if v[indx].is_fixed
                if haskey(indx_to_const_loc, indx)
                    const_loc = indx_to_const_loc[indx]
                    expr.nd[i] = NodeData(JuMP._Derivatives.VALUE, const_loc, nd.parent)
                else
                    push!(const_values, v[indx].lower_bound)
                    num_constants += 1
                    indx_to_const_loc[indx] = num_constants
                    expr.nd[i] = NodeData(nd.nodetype, num_constants, nd.parent)
                end
                f.isnumber[i] = true
            end
        end
    end

    return nothing
end

function eliminate_fixed_variables!(f::BufferedNonlinearFunction{N,T}, v::Vector{VariableInfo}) where {N,T}
    eliminate_fixed_variables!(f.expr, v)
end

function forward_pass!(x::Evaluator, d::NonlinearExpression{V}) where V
    for i = 1:d.dependent_subexpression_count
        !prior_eval(x, i) && forward_pass!(x, x.subexpressions[i])
    end
    load_subexprs!(d.b, x)
    fprop!(Relax, d.g, d.b)
    return
end

function forward_pass!(x::Evaluator, d::BufferedNonlinearFunction{V}) where V
    forward_pass!(x, d.expr)
    d.has_value = true
    d.last_past_reverse = false
    return
end

function reverse_pass!(evaluator::Evaluator, d::NonlinearExpression{V}) where V
    return rprop!(Relax, d.g, d.b)
end

function reverse_pass!(evaluator::Evaluator, d::CacheedNonlinearFunction{V}) where V
    d.last_past_reverse = true
    set_intersect_value!(d.expr, Interval(d.lower_bound, d.upper_bound))
    return reverse_pass!(evaluator, d.expr)
end
