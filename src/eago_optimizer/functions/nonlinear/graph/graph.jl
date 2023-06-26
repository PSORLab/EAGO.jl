# TODO: Each graph representation is assumed to be static... so

"""
$(TYPEDEF)

Abstract supertype for generic directed graph structure.
"""
abstract type AbstractDirectedGraph end
abstract type AbstractDirectedAcyclicGraph <: AbstractDirectedGraph end

const AbstractDG = AbstractDirectedGraph
const AbstractDAG = AbstractDirectedAcyclicGraph

@enum Linearity LIN_CONSTANT LIN_LINEAR LIN_PIECEWISE_LINEAR LIN_NONLINEAR
@enum VariableType VT_BIN VT_INT VT_CONT

function _variable_count(g::AbstractDG)::Int
    error("Variable count not defined for graph type = $(typeof(g))")
end

# Added id field to MOI OperatorRegistry
struct OperatorRegistry
    univariate_operators::Vector{Symbol}
    univariate_operator_id::Vector{Symbol}
    univariate_operator_to_id::Dict{Symbol,Int}
    univariate_user_operator_start::Int
    registered_univariate_operators::Vector{MOINL._UnivariateOperator}
    multivariate_operators::Vector{Symbol}
    multivariate_id::Vector{Symbol}
    multivariate_operator_to_id::Dict{Symbol,Int}
    multivariate_user_operator_start::Int
    registered_multivariate_operators::Vector{MOINL._MultivariateOperator}
    logic_operators::Vector{Symbol}
    logic_operator_id::Vector{Symbol}
    logic_operator_to_id::Dict{Symbol,Int}
    comparison_operators::Vector{Symbol}
    comparison_operator_id::Vector{Symbol}
    comparison_operator_to_id::Dict{Symbol,Int}
end
function OperatorRegistry()
    return OperatorRegistry(
        Symbol[],
        Symbol[],
        Dict{Symbol,Int}(),
        0,
        MOINL._UnivariateOperator[],
        Symbol[],
        Symbol[],
        Dict{Symbol,Int}(),
        0,
        MOINL._MultivariateOperator[],
        Symbol[],
        Symbol[],
        Dict{Symbol,Int}(),
        Symbol[],
        Symbol[],
        Dict{Symbol,Int}()
    )
end

include(joinpath(@__DIR__, "expressions.jl"))
include(joinpath(@__DIR__, "node.jl"))
include(joinpath(@__DIR__, "utilities.jl"))
include(joinpath(@__DIR__, "graphs", "directed_tree.jl"))