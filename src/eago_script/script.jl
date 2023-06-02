# Copyright (c) 2018: Matthew Wilhelm & Matthew Stuber.
# This code is licensed under MIT license (see LICENSE.md for full details)
#############################################################################
# EAGO
# A development environment for robust and global optimization
# See https://github.com/PSORLab/EAGO.jl
#############################################################################
# src/eago_script/script.jl
# A module used to manipulate script function inputs.
#############################################################################

module Script

    import MathOptInterface
    import MathOptInterface.Nonlinear: DEFAULT_MULTIVARIATE_OPERATORS, DEFAULT_UNIVARIATE_OPERATORS
    using MathOptInterface: AbstractOptimizer, features_available, initialize, NLPBlockData
    const MOINL = MathOptInterface.Nonlinear
    using SparseArrays: rowvals, nzrange, SparseMatrixCSC, spzeros
    import Base: afoldl, getindex, iterate
    import JuMP: NLPEvaluator
    import Cassette: @context, Cassette.overdub, Cassette.prehook
    import Cassette
    #import CodeTransformation: addmethod! TODO: add this later
    import ForwardDiff: derivative, gradient!
    export dag_flattening!

    include("codetransformation.jl") # TODO: delete this later
    #include("scrubber.jl")
    include("substitute.jl")
    include("patterns.jl")
    include("tracer.jl")
    include("loader.jl")
end
