# Copyright (c) 2018: Matthew Wilhelm & Matthew Stuber.
# This work is licensed under the Creative Commons Attribution-NonCommercial-
# ShareAlike 4.0 International License. To view a copy of this license, visit
# http://creativecommons.org/licenses/by-nc-sa/4.0/ or send a letter to Creative
# Commons, PO Box 1866, Mountain View, CA 94042, USA.
#############################################################################
# EAGO
# A development environment for robust and global optimization
# See https://github.com/PSORLab/EAGO.jl
#############################################################################
# src/EAGO.jl
# The main file for EAGO.
#############################################################################

__precompile__()

module EAGO

    import MathOptInterface

    using Reexport, Cassette, IntervalArithmetic, NumericIO, DocStringExtensions

    using JuMP
    import JuMP._Derivatives: operators, NodeData
    using JuMP._Derivatives: univariate_operators,
                             univariate_operator_to_id
    using Ipopt, GLPK

    using DataStructures: OrderedDict, BinaryMinMaxHeap, popmin!, popmax!, top
    using SparseArrays: SparseMatrixCSC, spzeros, rowvals, nzrange, nonzeros, sparse
    using LinearAlgebra: eigmin, norm

    import IntervalArithmetic: mid

    @reexport using McCormick
    @reexport using ReverseMcCormick

    const MOI = MathOptInterface

    const SV = MOI.SingleVariable
    const SAF = MOI.ScalarAffineFunction{Float64}
    const SQF = MOI.ScalarQuadraticFunction{Float64}
    const VECOFVAR = MOI.VectorOfVariables

    const SQT = MOI.ScalarQuadraticTerm{Float64}
    const SAT = MOI.ScalarAffineTerm{Float64}

    const LT = MOI.LessThan{Float64}
    const GT = MOI.GreaterThan{Float64}
    const ET = MOI.EqualTo{Float64}
    const IT = MOI.Interval{Float64}
    const ZO = MOI.ZeroOne
    const SOC = MOI.SecondOrderCone

    const VI = MOI.VariableIndex
    const CI = MOI.ConstraintIndex

    const SCoefC = MOI.ScalarCoefficientChange
    const SConsC = MOI.ScalarConstantChange
    const MOIU = MOI.Utilities

    const LT_ZERO = LT(0.0)

    # Add storage types for EAGO optimizers
    export NodeBB, get_history, get_lower_bound, get_upper_bound, get_lower_time,
           get_upper_time, get_preprocess_time, get_postprocess_time, get_lower_bound, get_solution_time,
           get_iteration_number, get_node_count, get_absolute_gap, get_relative_gap

    # map/reduce nonallocating no bounds checking map-reduce like utilities
    include("eago_optimizer/unsafe_utilities.jl")

    # creates a context that removes domain violations when constructing bounds
    include("eago_optimizer/guarded_context.jl")

    # defines structure used to store node in stack
    include("eago_optimizer/node_bb.jl")

    # load internal storage functions
    include("eago_optimizer/functions/functions.jl")

    #include("eago_optimizer/evaluator/evaluator.jl")

    # defines structure used to store information at each iteration of global optimize
    include("eago_optimizer/logging/log.jl")

    # defines the optimizer structures
    include("eago_optimizer/optimizer.jl")

    # defines routines to add variables and single variable constraints
    include("eago_optimizer/variables.jl")

    # defines routines to add saf, sqf, and nlp block constraints
    include("eago_optimizer/moi_constraints.jl")

    # functions which print information to console
    include("eago_optimizer/display.jl")

    #
    include("eago_optimizer/relax.jl")
    include("eago_optimizer/bound.jl")

    #
    include("eago_optimizer/domain_reduction.jl")

    #
    include("eago_optimizer/parse.jl")

    #
    include("eago_optimizer/logging/log_iteration.jl")

    #
    include("eago_optimizer/optimize/optimize.jl")

    # import the script solving utilities
    include("eago_script/script.jl")

    # routines for solving SIPs
    export SIP_Options, SIP_Result, explicit_sip_solve
    include("eago_semiinfinite/semi_infinite.jl")

end
