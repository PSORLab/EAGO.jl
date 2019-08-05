module Script

    using MathOptInterface: AbstractOptimizer
    import Base: afoldl, getindex
    import JuMP: _NLPData, _NonlinearExprData, _Derivatives.univariate_operator_to_id,
                 _Derivatives.operator_to_id, _Derivatives.comparison_operator_to_id,
                 _Derivatives.CALL, _Derivatives.CALLUNIVAR, _Derivatives.COMPARISON,
                 _Derivatives.VARIABLE, _Derivatives.VALUE, _Derivatives.SUBEXPRESSION,
                 _Derivatives.NodeType, _Derivatives.NodeData,
                 _Derivatives.USER_OPERATOR_ID_START, _Derivatives.USER_UNIVAR_OPERATOR_ID_START
    import Cassette: @context, Cassette.overdub
    #import CodeTransformation: addmethod! TODO: add this later
    import ForwardDiff: derivative, GradientConfig, gradient!

    include("codetransformation.jl") # TODO: delete this later
    include("scrubber.jl")
    include("tracer.jl")
    #include("loader.jl")
end
