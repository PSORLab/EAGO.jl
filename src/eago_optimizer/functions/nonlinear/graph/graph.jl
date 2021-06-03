#=
TODO: Each graph representation is assumed to be static... so
=#

abstract type AbstractDirectedGraph end
abstract type AbstractDirectedAcyclicGraph <: AbstractDirectedGraph end

const AbstractDG = AbstractDirectedGraph
const AbstractDAG = AbstractDirectedAcyclicGraph

@enum Linearity LIN_CONSTANT LIN_LINEAR LIN_PIECEWISE_LINEAR LIN_NONLINEAR
@enum VariableType VT_BIN VT_INT VT_CONT

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

include(joinpath(@__DIR__, "expressions.jl"))
include(joinpath(@__DIR__, "node.jl"))
include(joinpath(@__DIR__, "utilities.jl"))
include(joinpath(@__DIR__, "graphs", "directed_tree.jl"))