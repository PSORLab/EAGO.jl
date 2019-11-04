module Script

    using MathOptInterface: AbstractOptimizer, features_available, initialize, NLPBlockData
    using SparseArrays: rowvals, nzrange, SparseMatrixCSC, spzeros
    import Base: afoldl, getindex, iterate
    import JuMP: _NLPData, _NonlinearExprData
    import JuMP._Derivatives: univariate_operator_to_id, operator_to_id,
                              comparison_operator_to_id, CALL, CALLUNIVAR,
                              COMPARISON, VARIABLE, MOIVARIABLE, VALUE,
                              SUBEXPRESSION, PARAMETER, NodeType, NodeData,
                              USER_OPERATOR_ID_START, USER_UNIVAR_OPERATOR_ID_START,
                              UserOperatorRegistry, adjmat
    import Cassette: @context, Cassette.overdub, Cassette.prehook
    #import CodeTransformation: addmethod! TODO: add this later
    import ForwardDiff: derivative, gradient!
    export dag_flattening!

    include("codetransformation.jl") # TODO: delete this later
    include("scrubber.jl")
    include("substitute.jl")
    include("patterns.jl")
    include("tracer.jl")
    include("loader.jl")
end
