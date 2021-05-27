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
include(joinpath(@__DIR__, "constant", "constant.jl"))
include(joinpath(@__DIR__, "composite_relax", "composite_relax.jl"))

@enum(RelaxType, STD_RELAX, MC_AFF_RELAX, MC_ENUM_RELAX)

_set_has_value!(d, v) = nothing

"""
$(TYPEDEF)

Stores a general quadratic function with a buffer.
"""
mutable struct NonlinearExpression{V,S} <: AbstractEAGOConstraint
    g::DirectedTree{Float64}
    relax_cache::RelaxCache{V,S}
    has_value::Bool
    last_reverse::Bool
    lower_bound::S
    upper_bound::S
end
function NonlinearExpression()
    g = DirectedTree{Float64}()
    c = RelaxCache{MC{1,NS},Float64}()
    return NonlinearExpression{MC{1,NS},Float64}(g, c, false, false, -Inf, Inf)
end

function NonlinearExpression!(sub::Union{JuMP._SubexpressionStorage,JuMP._FunctionStorage},
                              b::MOI.NLPBoundsPair, sub_sparsity::Dict{Int,Vector{Int}},
                              subexpr_indx::Int,
                              subexpr_linearity::Vector{JuMP._Derivatives.Linearity},
                              op::OperatorRegistry, parameter_values,
                              tag::T; is_sub::Bool = false) where T
    g = DirectedTree{Float64}(sub, op, sub_sparsity, subexpr_linearity, parameter_values)
    grad_sparsity = _sparsity(g,1)
    n = length(grad_sparsity)
    if is_sub
        sub_sparsity[subexpr_indx] = copy(grad_sparsity) # updates subexpression sparsity dictionary
    end
    c = RelaxCache{MC{n,T},Float64}()
    initialize!(c, g)
    return NonlinearExpression{MC{n,T},Float64}(g, c, false, false, b.lower, b.upper)
end

@inline _has_value(d::NonlinearExpression) = d.has_value
@inline _dep_subexpr_count(d::NonlinearExpression) = _dep_subexpr_count(d.g)
@inline _set_has_value!(d::NonlinearExpression, v::Bool) = (d.has_value = v; return )
@inline _set_last_reverse!(d::NonlinearExpression, v::Bool) = (d.last_reverse = v; return )
@inline function _set_variable_storage!(d::NonlinearExpression, v::VariableValues{S}) where S<:Real
    d.relax_cache.v = v
    return
end
@inbounds _sparsity(d::NonlinearExpression) = _sparsity(d.g, 1)
@inbounds _set(d::NonlinearExpression{V,S}) where {V,S<:Real} = _set(d.relax_cache, 1)
@inbounds _num(d::NonlinearExpression{V,S}) where {V,S<:Real} = _num(d.relax_cache, 1)
@inbounds _is_num(d::NonlinearExpression) = _is_num(d.relax_cache, 1)

"""
$(TYPEDEF)

Stores a general nonlinear function with a buffer represented by the sum of a tape
and a scalar affine function.
"""
mutable struct BufferedNonlinearFunction{V,S} <: AbstractEAGOConstraint
    ex::NonlinearExpression{V,S}
    saf::SAF
end
function BufferedNonlinearFunction()
    ex = NonlinearExpression()
    saf = SAF(SAT[], 0.0)
    return BufferedNonlinearFunction{MC{1,NS},Float64}(ex, saf)
end

function BufferedNonlinearFunction(f::JuMP._FunctionStorage, b::MOI.NLPBoundsPair,
                                   sub_sparsity::Dict{Int,Vector{Int}},
                                   subexpr_lin::Vector{JuMP._Derivatives.Linearity},
                                   op::OperatorRegistry, parameter_values,
                                   tag::T) where T <: RelaxTag

    ex = NonlinearExpression!(f, b, sub_sparsity, -1, subexpr_lin, op, parameter_values, tag)
    n = length(_sparsity(ex.g, 1))
    saf = SAF(SAT[SAT(0.0, VI(i)) for i = 1:n], 0.0)
    return BufferedNonlinearFunction{MC{n,T},Float64}(ex, saf)
end

@inline _set_last_reverse!(d::BufferedNonlinearFunction, v::Bool) = _set_last_reverse!(d.ex, v)
@inline function _set_variable_storage!(d::BufferedNonlinearFunction, v::VariableValues{S}) where S<:Real
    _set_variable_storage!(d.ex, v)
end

@inline _has_value(d::BufferedNonlinearFunction) where {V,S<:Real} = _has_value(d.ex)
@inline _dep_subexpr_count(d::BufferedNonlinearFunction) = _dep_subexpr_count(d.ex)
@inline _set_has_value!(d::BufferedNonlinearFunction, v::Bool) = _set_has_value!(d.ex, v)
@inline _sparsity(d::BufferedNonlinearFunction) = _sparsity(d.ex)
@inline _set(d::BufferedNonlinearFunction{V,S}) where {V,S<:Real} = _set(d.ex)
@inline _num(d::BufferedNonlinearFunction{V,S}) where {V,S<:Real} = _num(d.ex)
@inline _lower_bound(d::BufferedNonlinearFunction{V,S}) where {V,S<:Real} = d.ex.lower_bound
@inline _upper_bound(d::BufferedNonlinearFunction{V,S}) where {V,S<:Real} = d.ex.upper_bound
# returns the interval bounds associated with the set
@inline _interval(d::BufferedNonlinearFunction{V,S}) where {V,S<:Real} = Interval{S}(_set(d))
@inline _is_num(d::BufferedNonlinearFunction) = _is_num(d.ex)


"""
    Evaluator

MOI.AbstractNLPEvaluator for calculating relaxations of nonlinear terms.

Checks that the resulting value should be a number...

$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct Evaluator <: MOI.AbstractNLPEvaluator
    user_operators::OperatorRegistry = OperatorRegistry()
    has_user_mv_operator::Bool = false
    num_mv_buffer::Vector{Float64} = Float64[]
    parameter_values::Vector{Float64} = Float64[]
    node::NodeBB = NodeBB()
    variable_values::VariableValues{Float64} = VariableValues{Float64}()
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
    relax_type::RelaxType                       = STD_RELAX
end

"""
$(FUNCTIONNAME)

Sets the current node in the Evaluator structure.
"""
function set_node!(d::Evaluator, n::NodeBB)
    d.node = NodeBB(n)
    @inbounds for i = 1:length(n)
        vi = d.variable_values.node_to_variable_map[i]
        d.variable_values.lower_variable_bounds[vi] = n.lower_variable_bounds[i]
        d.variable_values.upper_variable_bounds[vi] = n.upper_variable_bounds[i]
    end
    fill!(d.subexpressions_eval, false)
    d.is_first_eval = true
    return nothing
end

function retrieve_node(d::Evaluator)
    n = d.node
    nv_map = d.variable_values.node_to_variable_map
    return NodeBB(copy(d.variable_values.lower_variable_bounds[nv_map]),
                  copy(d.variable_values.upper_variable_bounds[nv_map]),
                  copy(n.is_integer), n.continuous,
                  n.lower_bound, n.upper_bound, n.depth, n.id,
                  n.branch_direction, n.last_branch, n.branch_extent)
end
@inline function _get_x!(::Type{BranchVar}, out::Vector{Float64}, d::Evaluator)
    return _get_x!(BranchVar, out, d.variable_values)
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

function eliminate_fixed_variables!(f::BufferedNonlinearFunction{V,S}, v::Vector{VariableInfo}) where {V,S}
    eliminate_fixed_variables!(f.ex, v)
end

function forward_pass!(x::Evaluator, d::NonlinearExpression{V}) where V  # Hold reference to subexpressions in DAG?
    # Fix subexpression code...
    #for i = 1:_dep_subexpr_count(d)
    #    !prior_eval(x, i) && forward_pass!(x, x.subexpressions[i])
    #end
    #_load_subexprs!(d.relax_cache, x.subexpressions)
    if x.relax_type == STD_RELAX
        fprop!(Relax(), d.g, d.relax_cache)
    elseif x.relax_type == MC_AFF_RELAX
        fprop!(RelaxAA(), d.g, d.relax_cache)
    elseif x.relax_type == MC_ENUM_RELAX
        fprop!(RelaxMulEnum(), d.g, d.relax_cache)
    end
    return
end

function forward_pass!(x::Evaluator, d::BufferedNonlinearFunction{V,S}) where {V,S}
    forward_pass!(x, d.ex)
    _set_has_value!(d, true)
    _set_last_reverse!(d, false)
    return
end

function rprop!(::Relax, x::Evaluator, d::NonlinearExpression{V}) where V
    return rprop!(Relax(), d.g, d.relax_cache)
end
function rprop!(::Relax, x::Evaluator, d::BufferedNonlinearFunction{V,S}) where {V,S}
    _set_last_reverse!(d, true)
    return rprop!(Relax(), x, d.ex)
end
