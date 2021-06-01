#=
TODO: Each graph representation is assumed to be static... so
=#

abstract type AbstractDirectedGraph end
abstract type AbstractDirectedAcyclicGraph <: AbstractDirectedGraph end

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

"""
    _rev_sparsity

Return the index of the ith variable at node k.
"""
@inline _rev_sparsity(g::AbstractDG, i::Int, k::Int) = error("_rev_sparsity not defined for g::$(typeof(g))")

# added id field from JuMP UserOperatorRegistry, expect more extensive changes in future.
struct OperatorRegistry
    multivariate_id::Vector{Symbol}
    multivariate_operator_to_id::Dict{Symbol,Int}
    multivariate_operator_evaluator::Vector{MOI.AbstractNLPEvaluator}
    univariate_operator_id::Vector{Symbol}
    univariate_operator_to_id::Dict{Symbol,Int}
    univariate_operator_f::Vector{Any}
    univariate_operator_fprime::Vector{Any}
    univariate_operator_fprimeprime::Vector{Any}
end
function OperatorRegistry()
    return OperatorRegistry(
        Symbol[],
        Dict{Symbol,Int}(),
        MOI.AbstractNLPEvaluator[],
        Symbol[],
        Dict{Symbol,Int}(),
        [],
        [],
        [],
    )
end

include(joinpath(@__DIR__, "utilities.jl"))
include(joinpath(@__DIR__, "expressions.jl"))
include(joinpath(@__DIR__, "abstract_node.jl"))
include(joinpath(@__DIR__, "abstract_cache.jl"))
include(joinpath(@__DIR__, "graphs", "directed_tree.jl"))

const ALLGRAPHS = Union{DAT}

# node property access functions that can be defined at abstract type
@inline _node(g::ALLGRAPHS, i)           = @inbounds getindex(_nodes(g), i)
@inline _variable(g::ALLGRAPHS, i)       = @inbounds getindex(_variables(g), i)
@inline _variable_type(g::ALLGRAPHS, i)  = @inbounds getindex(_variable_types(g), i)
@inline function _constant_value(g::ALLGRAPHS, i)
    @inbounds getindex(_constant_values(g), i)
end
@inline function _parameter_value(g::ALLGRAPHS, i)
    @inbounds getindex(_parameter_values(g), i)
end

@inline _node_class(g::ALLGRAPHS, i)      = _node_class(_node(g, i))
@inline _ex_type(g::ALLGRAPHS, i)         = _ex_type(_node(g, i))
@inline _first_index(g::ALLGRAPHS, i)     = _first_index(_node(g, i))
@inline _secondary_index(g::ALLGRAPHS, i) = _secondary_index(_node(g, i))
@inline _arity(g::ALLGRAPHS, i)           = _arity(_node(g, i))
@inline _children(g::ALLGRAPHS, i)        = _children(_node(g, i))
@inline _child(g::ALLGRAPHS, i, j)        = _child(_node(g, j), i)

"""
    _node_count

Number of nodes in graph g.
"""
@inline _node_count(g::ALLGRAPHS) = length(_nodes(g))

"""
    _variable_count

Number of variables in graph g.
"""
@inline _variable_count(g::ALLGRAPHS) = length(_variable(g))

"""
    _constant_count

Number of constants in graph g.
"""
@inline _constant_count(g::ALLGRAPHS) = length(_constant_count(g))
