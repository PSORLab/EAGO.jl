# Copyright (c) 2018: Matthew Wilhelm & Matthew Stuber.
# This code is licensed under MIT license (see LICENSE.md for full details)
#############################################################################
# Disciplined.jl
# DAG Representation & Utilities for Optimization
# See https://github.com/PSORLab/EAGO.jl
#############################################################################
# src/functions/nonlinear/g/c.jl
# Defines the AbstractCache and AbstractCacheAttribute types along with
# fallback functions for undefined.
#############################################################################

"""
    AbstractCache

Abstract supertype used for information storage object the directed acyclic graph.
"""
abstract type AbstractCache end

"""
    AbstractCacheAttribute

Abstract supertype used for attributes stored in a cache.
"""
abstract type AbstractCacheAttribute end

"""
    initialize!

Initializes the cache `c` to correspond to directed graph `g`.
"""
function initialize!(c::AbstractCache, g::AbstractDG)
    error("No function initialize!(g,c) defined for g = $(typeof(g)) and c = $(typeof(c)).")
end

"""
    fprop!

Performs a forward walk on the `g::DirectedAcyclicGraph` to calculate the
attribute `t::AbstractCacheAttribute` stored in `c::AbstractCache`. An
optional fourth parameter `i::Int` indicates the node in the graph to evaluate
(evaluating all dependent nodes as necessary).
"""
function fprop!(::AbstractCacheAttribute, g::AbstractDG, c::AbstractCache)
    error("No function fprop!(t, g, c) defined for t = $(typeof(t)), g = $(typeof(g)) and c = $(typeof(c)).")
end
function fprop!(::AbstractCacheAttribute, g::AbstractDG, c::AbstractCache, i::Int)
    error("No function fprop!(t, g, c, i) defined for t = $(typeof(t)), g = $(typeof(g)) and c = $(typeof(c)).")
end

"""
    rprop!

Performs a reverse walk on the `g::DirectedAcyclicGraph` to calculate the
attribute `t::AbstractCacheAttribute` stored in `c::AbstractCache`. An
optional fourth parameter `i::Int` indicates the node in the graph to begin the
`reverse` evaluation from.
"""
function rprop!(t::AbstractCacheAttribute, g::AbstractDG, c::AbstractCache)
    error("No function rprop!(t, g, b) defined for t::AbstractCacheAttribute
           = $(typeof(t)), g::AbstractDG = $(typeof(g)) and c::AbstractCache = $(typeof(c)).")
end
function rprop!(t::AbstractCacheAttribute, g::AbstractDG, c::AbstractCache, i::Int)
    error("No function rprop!(t, g, b) defined for t::AbstractCacheAttribute
           = $(typeof(t)), g::AbstractDG = $(typeof(g)) and c::AbstractCache = $(typeof(c)).")
end

"""
    _is_discovered(c, i)

Check if a node was discovered by search.
"""
function _is_discovered(t::AbstractCacheAttribute, c::AbstractCache, i::Int)
    error("No function _is_discovered(t, c, i) defined for t::AbstractCacheAttribute = $t and c::AbstractCache = $c.")
end

"""
    _is_locked(t::AbstractCacheAttribute, c::AbstractCache, i::Int)

Defines whether the attribute associated with each node is locked.
A lock node has been discovered.  A discovered not has not necessarily been locked.
"""
function _is_locked(t::AbstractCacheAttribute, c::AbstractCache, i::Int)
    error("No function _is_locked(t, c, i) defined for t = $t and c = $c.")
end

"""
    discover!

Specifies that a node has been discovered by a search.
"""
function discover!(t::AbstractCacheAttribute, c::AbstractCache, i::Int)
    error("No function discovered!(t, c, i) defined for t = $t and c = $c.")
end

_discoverable(attr, g, c, i) = !_is_discovered(attr, g, c, i)
_unlocked(attr, g, c, i) = !_is_locked(attr, g, c, i)


Base.@kwdef mutable struct VariableValues{T<:Real}
    x::Vector{T}                      = T[]
    lower_variable_bounds::Vector{T}  = T[]
    upper_variable_bounds::Vector{T}  = T[]
    node_to_variable_map::Vector{Int} = Int[]
    variable_to_node_map::Vector{Int} = Int[]
end

@inline _val(b::VariableValues{T}, i::Int) where T = @inbounds b.x[i]
@inline _lbd(b::VariableValues{T}, i::Int) where T = @inbounds b.lower_variable_bounds[i]
@inline _ubd(b::VariableValues{T}, i::Int) where T = @inbounds b.upper_variable_bounds[i]

# Define standard forward and reverse propagation to switch of expression definitions
# for expressions.
function binary_switch(ids, exprs; is_rev = true)
    if length(exprs) <= 3
        if is_rev
            out = Expr(:if, Expr(:call, :(==), :id, ids[1]),
                       :(fprop!(T, $(exprs[1]), g, b, k)))
        else
            out = Expr(:if, Expr(:call, :(==), :id, ids[1]),
                       :(rprop!(T, $(exprs[1]), g, b, k)))
        end
        if length(exprs) > 1
            push!(out.args, binary_switch(ids[2:end], exprs[2:end]))
        end
        return out
    else
        mid = length(exprs) >>> 1
        return Expr(:if, Expr(:call, :(<=), :id, ids[mid]),
                         binary_switch(ids[1:mid], exprs[1:mid]),
                         binary_switch(ids[mid+1:end], exprs[mid+1:end]))
    end
end

const ALL_ATM_EVAL = [ALL_ATOM_DICT[i] for i in ALL_ATOM_TYPES]
f_switch = binary_switch(ALL_ATOM_TYPES, ALL_ATM_EVAL, is_rev = false)
r_switch = binary_switch(ALL_ATOM_TYPES, ALL_ATM_EVAL, is_rev = true)
@eval @inline function fprop!(::Type{T}, ::Type{Expression}, g::AbstractDG, c::AbstractCache , k::Int) where T<:AbstractCacheAttribute
    id = _ex_type(g, k)
    $f_switch
end
@eval @inline function rprop!(::Type{T}, ::Type{Expression}, g::AbstractDG, c::AbstractCache , k::Int) where T<:AbstractCacheAttribute
    id = _ex_type(g, k)
    $r_switch
end
