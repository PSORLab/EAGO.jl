#=
TODO: Each graph representation is assumed to be static... so
=#

abstract type AbstractDirectedGraph{S<:Real} end
abstract type AbstractDirectedAcyclicGraph{S<:Real} <: AbstractDirectedGraph{S} end

const AbstractDG = AbstractDirectedGraph
const AbstractDAG = AbstractDirectedAcyclicGraph

@enum Linearity LIN_CONSTANT LIN_LINEAR LIN_PIECEWISE_LINEAR LIN_NONLINEAR
@enum VariableType VT_BIN VT_INT VT_CONT

"""
    _nodes

Access node list for graph g.
"""
_nodes(g::AbstractDG) = error("_nodes(g) not defined for g::$(typeof(g)).")

"""
    _variable

Access variable list for graph g.
"""
_variable(g::AbstractDG) = error("_variables(g) not defined for g::$(typeof(g)).")

"""
    _variable_types

Access variable type list for graph g.
"""
_variable_types(g::AbstractDG, i) = error("_variable_types(g) not defined for g::$(typeof(g)).")

"""
    _constant_values

Access constant values list for graph g.
"""
_constant_values(g::AbstractDG, i) = error("_constant_values(g) not defined for g::$(typeof(g)).")

"""
    _user_univariate_operator

Retreive the ith univariate user function evaluator stored in the graph.
"""
_user_univariate_operator(g::AbstractDG, i) = error("_user_univariate_operator(g,i) not defined for g::$(typeof(g)).")

"""
    _user_univariate_operator

Retreive the ith multivariate user function evaluator stored in the graph.
"""
_user_multivariate_operator(g::AbstractDG, i) = error("_user_multivariate_operator(g,i) not defined for g::$(typeof(g)).")

# node property access functions that can be defined at abstract type
@inline _node(g::AbstractDG, i)           = @inbounds getindex(_nodes(g), i)
@inline _variable(g::AbstractDG, i)       = @inbounds getindex(_variables(g), i)
@inline _variable_type(g::AbstractDG, i)  = @inbounds getindex(_variable_types(g), i)
@inline function _constant_value(g::AbstractDG{S}, i) where S <: Real
    @inbounds getindex(_constant_values(g), i)
end

@inline _node_type(g::AbstractDG, i)       = _node_type(_node(g, i))
@inline _expr_type(g::AbstractDG, i)       = _expr_type(_node(g, i))
@inline _index(g::AbstractDG, i)           = _index(_node(g, i))
@inline _secondary_index(g::AbstractDG, i) = _secondary_index(_node(g, i))
@inline _arity(g::AbstractDG, i)           = _arity(_node(g, i))
@inline _children(g::AbstractDG, i)        = _children(_node(g, i))
@inline _child(g::AbstractDG, i, j)        = _child(_node(g, j), i)

"""
    _node_count

Number of nodes in graph g.
"""
@inline _node_count(g::AbstractDG) = length(_nodes(g))

"""
    _variable_count

Number of variables in graph g.
"""
@inline _variable_count(g::AbstractDG) = length(_variable(g))

"""
    _constant_count

Number of constants in graph g.
"""
@inline _constant_count(g::AbstractDG) = length(_constant_count(g))

include(joinpath(@__DIR__, "utilities.jl"))
include(joinpath(@__DIR__, "expressions.jl"))
include(joinpath(@__DIR__, "abstract_node.jl"))
include(joinpath(@__DIR__, "abstract_cache.jl"))
include(joinpath(@__DIR__, "graphs", "directed_tree.jl"))
