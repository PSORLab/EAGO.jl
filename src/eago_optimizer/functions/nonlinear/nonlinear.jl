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

include(joinpath(@__DIR__, "register_special.jl"))
include(joinpath(@__DIR__, "graph", "graph.jl"))
include(joinpath(@__DIR__, "constant", "constant.jl"))
include(joinpath(@__DIR__, "composite_relax", "composite_relax.jl"))
include(joinpath(@__DIR__, "apriori_relax", "apriori_relax.jl"))

@enum(RelaxType, STD_RELAX, MC_AFF_RELAX, MC_ENUM_RELAX)

_set_has_value!(d, v) = nothing

"""
$(TYPEDEF)

Stores a general quadratic function with a buffer.
"""
mutable struct NonlinearExpression{V,N,T<:RelaxTag} <: AbstractEAGOConstraint
    g::DirectedTree
    relax_cache::RelaxCache{V,N,T}
    has_value::Bool
    last_reverse::Bool
    lower_bound::Float64
    upper_bound::Float64
    grad_sparsity::Vector{Int}
end
function NonlinearExpression()
    g = DirectedTree()
    c = RelaxCache{MC{1,NS},1,NS}()
    return NonlinearExpression{MC{1,NS},1,NS}(g, c, false, false, -Inf, Inf, Int[])
end

relax_info(s::Relax, n::Int, t::T) where T = MC{n,T}
function NonlinearExpression!(rtype::S, sub::Union{JuMP._SubexpressionStorage,JuMP._FunctionStorage},
                              b::MOI.NLPBoundsPair, sub_sparsity::Dict{Int,Vector{Int}},
                              subexpr_indx::Int,
                              subexpr_linearity::Vector{JuMP._Derivatives.Linearity},
                              op::OperatorRegistry, parameter_values,
                              tag::T, use_apriori_flag::Bool; is_sub::Bool = false) where {S,T}
    g = DirectedTree(sub, op, sub_sparsity, subexpr_linearity, parameter_values)
    grad_sparsity = _sparsity(g, 1)
    n = length(grad_sparsity)
    if is_sub
        sub_sparsity[subexpr_indx] = copy(grad_sparsity) # updates subexpression sparsity dictionary
    end
    V = relax_info(rtype, n, tag)
    c = RelaxCache{V,n,T}()
    c.use_apriori_mul = use_apriori_flag
    initialize!(c, g)
    return NonlinearExpression{V,n,T}(g, c, false, false, b.lower, b.upper, grad_sparsity)
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
@inbounds _set(d::NonlinearExpression{V}) where V = _set(d.relax_cache, 1)
@inbounds _num(d::NonlinearExpression{V}) where V = _num(d.relax_cache, 1)
@inbounds _is_num(d::NonlinearExpression) = _is_num(d.relax_cache, 1)

"""
$(TYPEDEF)

Stores a general nonlinear function with a buffer represented by the sum of a tape
and a scalar affine function.
"""
mutable struct BufferedNonlinearFunction{V,N,T<:RelaxTag} <: AbstractEAGOConstraint
    ex::NonlinearExpression{V,N,T}
    saf::SAF
end
function BufferedNonlinearFunction()
    ex = NonlinearExpression()
    saf = SAF(SAT[], 0.0)
    return BufferedNonlinearFunction{MC{1,NS},1,NS}(ex, saf)
end

function BufferedNonlinearFunction(rtype::RELAX_ATTRIBUTE, f::JuMP._FunctionStorage, b::MOI.NLPBoundsPair,
                                   sub_sparsity::Dict{Int,Vector{Int}},
                                   subexpr_lin::Vector{JuMP._Derivatives.Linearity},
                                   op::OperatorRegistry, parameter_values,
                                   tag::T, use_apriori_flag::Bool) where T <: RelaxTag

    ex = NonlinearExpression!(rtype, f, b, sub_sparsity, -1, subexpr_lin, op, parameter_values, tag, use_apriori_flag)
    n = length(_sparsity(ex.g, 1))
    saf = SAF(SAT[SAT(0.0, VI(i)) for i = 1:n], 0.0)
    V = relax_info(rtype, n, tag)
    return BufferedNonlinearFunction{V,n,T}(ex, saf)
end

@inline _set_last_reverse!(d::BufferedNonlinearFunction{V,N,T}, v::Bool) where {V,N,T<:RelaxTag} = _set_last_reverse!(d.ex, v)
function _set_variable_storage!(d::BufferedNonlinearFunction{V,N,T}, v::VariableValues{Float64}) where {V,N,T<:RelaxTag}
    _set_variable_storage!(d.ex, v)
end

_has_value(d::BufferedNonlinearFunction) = _has_value(d.ex)
_dep_subexpr_count(d::BufferedNonlinearFunction) = _dep_subexpr_count(d.ex)
_set_has_value!(d::BufferedNonlinearFunction, v::Bool) = _set_has_value!(d.ex, v)
_sparsity(d::BufferedNonlinearFunction) = _sparsity(d.ex)
_set(d::BufferedNonlinearFunction{V,N,T}) where {V,N,T<:RelaxTag} = _set(d.ex)
_num(d::BufferedNonlinearFunction{V,N,T}) where {V,N,T<:RelaxTag} = _num(d.ex)
_lower_bound(d::BufferedNonlinearFunction{V,N,T}) where {V,N,T<:RelaxTag} = d.ex.lower_bound
_upper_bound(d::BufferedNonlinearFunction{V,N,T}) where {V,N,T<:RelaxTag} = d.ex.upper_bound
# returns the interval bounds associated with the set
_interval(d::BufferedNonlinearFunction{V,N,T}) where {V,N,T<:RelaxTag} = Interval{Float64}(_set(d))
_is_num(d::BufferedNonlinearFunction) = _is_num(d.ex)


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
    for i = 1:length(n)
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
                  n.lower_bound, n.upper_bound, n.depth, n.cont_depth, n.id,
                  n.branch_direction, n.last_branch, n.branch_extent)
end
@inline function _get_x!(::Type{BranchVar}, out::Vector{Float64}, d::Evaluator)
    return _get_x!(BranchVar, out, d.variable_values)
end
prior_eval(d::Evaluator, i::Int) = d.subexpressions_eval[i]

#=
Assumes the sparsities are sorted...
=#
function copy_subexpression_value!(k::Int, op::Int, subexpression::NonlinearExpression{V,MC{N1,T}},
                                   numvalued::Vector{Bool}, numberstorage::Vector{S}, setstorage::Vector{MC{N2,T}},
                                   cv_buffer::Vector{S}, cc_buffer::Vector{S},
                                   func_sparsity::Vector{Int}) where {V, N1, N2, S, T <: RelaxTag}

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

function eliminate_fixed_variables!(f::NonlinearExpression{V,N,T}, v::Vector{VariableInfo}) where {V,N,T<:RelaxTag}
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

function eliminate_fixed_variables!(f::BufferedNonlinearFunction{N,T}, v::Vector{VariableInfo}) where {N,T<:RelaxTag}
    eliminate_fixed_variables!(f.ex, v)
end

function f_init_prop!(t, g::DAT, c::RelaxCache, flag::Bool)
    if flag
        return f_init!(t, g, c)
    end
    return fprop!(t, g, c)
end
function forward_pass!(x::Evaluator, d::NonlinearExpression{V,N,T}) where {V,N,T<:RelaxTag}
    # Fix subexpression code...
    #for i = 1:_dep_subexpr_count(d)
    #    !prior_eval(x, i) && forward_pass!(x, x.subexpressions[i])
    #end
    #_load_subexprs!(d.relax_cache, x.subexpressions)
    (x.relax_type == STD_RELAX)    && (return f_init_prop!(Relax(), d.g, d.relax_cache, x.is_first_eval))
    (x.relax_type == MC_AFF_RELAX) && (return f_init_prop!(RelaxAA(d.grad_sparsity), d.g, d.relax_cache, x.is_first_eval))
    return f_init_prop!(RelaxMulEnum(d.grad_sparsity), d.g, d.relax_cache, x.is_first_eval)
end

function forward_pass!(x::Evaluator, d::BufferedNonlinearFunction{V,N,T}) where {V,N,T<:RelaxTag}
    forward_pass!(x, d.ex)
    _set_has_value!(d, true)
    _set_last_reverse!(d, false)
    return
end

rprop!(::Relax, x::Evaluator, d::NonlinearExpression{V,N,T}) where {V,N,T<:RelaxTag} = rprop!(Relax(), d.g, d.relax_cache)
function rprop!(::Relax, x::Evaluator, d::BufferedNonlinearFunction{V,N,T}) where {V,N,T<:RelaxTag}
    _set_last_reverse!(d, true)
    return rprop!(Relax(), x, d.ex)
end
