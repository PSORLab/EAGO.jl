__precompile__()

module EAGO

    import MathOptInterface

    using Reexport, Cassette, IntervalArithmetic, NumericIO, DocStringExtensions

    using JuMP
    import JuMP._Derivatives: operators, NodeData
    using JuMP._Derivatives: univariate_operators,
                             univariate_operator_to_id
    using Ipopt, GLPK

    using DataStructures: BinaryMinMaxHeap, popmin!, popmax!, top
    using SparseArrays: SparseMatrixCSC, spzeros, rowvals, nzrange, nonzeros, sparse
    using LinearAlgebra: eigmin, norm

    import IntervalArithmetic: mid

    @reexport using McCormick
    @reexport using ReverseMcCormick

    import Base: ImmutableDict, isless, isempty, eltype, copy, length

    const MOI = MathOptInterface
    const SAF = MOI.ScalarAffineFunction{Float64}
    const SAT = MOI.ScalarAffineTerm{Float64}
    const SQF = MOI.ScalarQuadraticFunction{Float64}
    const SQT = MOI.ScalarQuadraticTerm{Float64}
    const SV = MOI.SingleVariable
    const VECOFVAR = MOI.VectorOfVariables
    const LT = MOI.LessThan{Float64}
    const GT = MOI.GreaterThan{Float64}
    const ET = MOI.EqualTo{Float64}
    const IT = MOI.Interval{Float64}
    const ZO = MOI.ZeroOne
    const VI = MOI.VariableIndex
    const CI = MOI.ConstraintIndex
    const SCoefC = MOI.ScalarCoefficientChange
    const SConsC = MOI.ScalarConstantChange
    const MOIU = MOI.Utilities

    # Add storage types for EAGO optimizers
    export NodeBB, get_history, get_lower_bound, get_upper_bound, get_lower_time,
           get_upper_time, get_preprocess_time, get_postprocess_time, get_lower_bound, get_solution_time,
           get_iteration_number, get_node_count, get_absolute_gap, get_relative_gap

    include("eago_optimizer/unsafe_utilities.jl")
    include("eago_optimizer/guarded_context.jl")
    include("eago_optimizer/node_bb.jl")
    include("eago_optimizer/evaluator/evaluator.jl")
    include("eago_optimizer/logging.jl")
    include("eago_optimizer/optimizer.jl")
    include("eago_optimizer/variables.jl")
    include("eago_optimizer/constraints.jl")
    include("eago_optimizer/display.jl")
    include("eago_optimizer/relax.jl")
    include("eago_optimizer/domain_reduction.jl")
    include("eago_optimizer/subroutines.jl")
    include("eago_optimizer/optimize.jl")

    # Import the script solving utilities
    include("eago_script/script.jl")

    # Routines for solving SIPs
    export SIP_Options, SIP_Result, explicit_sip_solve
    include("eago_semiinfinite/semi_infinite.jl")

end
