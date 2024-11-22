# Copyright (c) 2018: Matthew Wilhelm, Robert Gottlieb, Dimitri Alston,
# Matthew Stuber, and the University of Connecticut (UConn).
# This code is licensed under the MIT license (see LICENSE.md for full details).
################################################################################
# EAGO
# A development environment for robust and global optimization.
# https://github.com/PSORLab/EAGO.jl
################################################################################
# src/EAGO.jl
# The main file for EAGO.
################################################################################

__precompile__()

module EAGO

    using MathOptInterface
    import MathOptInterface
    import MathOptInterface.Nonlinear: DEFAULT_MULTIVARIATE_OPERATORS, DEFAULT_UNIVARIATE_OPERATORS
    using Reexport, Requires, Cassette, IntervalArithmetic, DocStringExtensions,
          FastRounding, SpecialFunctions, Ipopt, Cbc, Printf, PrettyTables

    using JuMP
    import JuMP

    using DataStructures: OrderedDict, BinaryMinMaxHeap, popmin!, popmax!, top
    using SparseArrays: SparseMatrixCSC, spzeros, rowvals, nzrange, nonzeros, sparse, findnz
    using LinearAlgebra: eigmin, norm
    using Base: @propagate_inbounds
    import Base: setindex!, + , *, -, ^, /, zero, one, inv, log, log10, exp, exp10, isempty, min, max

    @reexport using McCormick
    import McCormick: relu, leaky_relu, maxsig, maxtanh, softplus, pentanh, sigmoid, bisigmoid,
                      softsign, gelu, swish, xabsx, logcosh, xlogx, erf, erfinv, erfc, erfcinv, 
                      xexpax, arh, param_relu, elu, selu, lower_bnd, upper_bnd, bnd, positive, negative

    @reexport using SpecialFunctions
    #@reexport using ReverseMcCormick

    using IntervalContractors
    #=
    using IntervalContractors: plus_rev, minus_rev, inv_rev,
       mul_rev, div_rev, power_rev,
       max_rev, min_rev,
        sqr_rev, sqrt_rev, abs_rev,
        exp_rev, exp2_rev, exp10_rev, expm1_rev,
        log_rev, log2_rev, log10_rev, log1p_rev,
        sin_rev, cos_rev, tan_rev,
        asin_rev, acos_rev, atan_rev,
        sinh_rev, cosh_rev, tanh_rev,
        asinh_rev, acosh_rev, atanh_rev
        =#

    const MOI = MathOptInterface
    const MOIU = MOI.Utilities
    const MOIB = MOI.Bridges
    const MOINL = MOI.Nonlinear
    const MOIRAD = MOINL.ReverseAD

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
    const INT = MOI.Integer
    const SOC = MOI.SecondOrderCone

    const VI = MOI.VariableIndex
    const CI = MOI.ConstraintIndex

    const SCoefC = MOI.ScalarCoefficientChange
    const SConsC = MOI.ScalarConstantChange

    const LT_ZERO = LT(0.0)

    # Add storage types for EAGO optimizers
    export NodeBB, get_history, get_lower_bound, get_upper_bound, get_lower_time,
           get_upper_time, get_preprocess_time, get_postprocess_time, get_lower_bound, get_solution_time,
           get_iteration_number, get_node_count, get_absolute_gap, get_relative_gap, SubSolvers, 
           EAGOModel, AuxiliaryVariableRef
    export auxiliary_variable, @auxiliary_variable, add_mimo_expression, @add_mimo_expression 

    export relu, leaky_relu, maxsig, maxtanh, softplus, pentanh, sigmoid, bisigmoid,
           softsign, gelu, swish, xabsx, logcosh, xlogx, erf, erfinv, erfc, erfcinv, 
           xexpax, arh, param_relu, elu, selu, lower_bnd, upper_bnd, bnd, positive, negative

    # Map/reduce nonallocating no bounds checking map-reduce like utilities
    include(joinpath(@__DIR__, "eago_optimizer", "debug_tools.jl"))
    include(joinpath(@__DIR__, "eago_optimizer", "utilities.jl"))

    # Creates a context that removes domain violations when constructing bounds
    #include("eago_optimizer/guarded_context.jl")

    include(joinpath(@__DIR__, "eago_optimizer", "types", "log.jl"))
    include(joinpath(@__DIR__, "eago_optimizer", "types", "variable_info.jl"))
    include(joinpath(@__DIR__, "eago_optimizer", "types", "node_bb.jl"))
    include(joinpath(@__DIR__, "eago_optimizer", "types", "extension.jl"))
    include(joinpath(@__DIR__, "eago_optimizer", "types", "incremental.jl"))
    include(joinpath(@__DIR__, "eago_optimizer", "types", "subsolver_block.jl"))

    # Loads internal storage functions
    include(joinpath(@__DIR__, "eago_optimizer", "functions", "functions.jl"))

    include(joinpath(@__DIR__, "eago_optimizer", "types", "global_optimizer.jl"))

    # Defines the optimizer structures
    include(joinpath(@__DIR__, "eago_optimizer", "optimizer.jl"))

    # Defines routines to add variable, saf, sqf, and nlp block constraints
    include(joinpath(@__DIR__, "eago_optimizer", "moi_wrapper.jl"))

    #
    include(joinpath(@__DIR__, "eago_optimizer", "optimize", "nonconvex", "relax.jl"))
    include(joinpath(@__DIR__, "eago_optimizer", "optimize", "nonconvex", "bound.jl"))

    #
    include(joinpath(@__DIR__, "eago_optimizer", "domain_reduction.jl"))

    #
    include(joinpath(@__DIR__, "eago_optimizer", "parse.jl"))

    #
    include(joinpath(@__DIR__, "eago_optimizer", "optimize", "optimize.jl"))

    # Imports the script solving utilities
    include(joinpath(@__DIR__, "eago_script", "script.jl"))

    # Routines for solving SIPs
    export SIPResult, SIPProblem, SIPCallback, SIPSubResult,
           sip_solve, SIPRes, SIPResRev, SIPHybrid,
           build_model, set_tolerance_inner!, set_tolerance!, get_disc_set,
           sip_llp!, sip_bnd!, sip_res!, get_sip_optimizer, check_convergence,
           LowerLevel1, LowerLevel2, LowerLevel3, LowerProblem, UpperProblem,
           ResProblem, AbstractSIPAlgo, AbstractSubproblemType
    include(joinpath(@__DIR__, "eago_semiinfinite", "semiinfinite.jl"))

    include(joinpath(@__DIR__, "subsolvers", "cbc.jl"))
    include(joinpath(@__DIR__, "subsolvers", "ipopt.jl"))
    function __init__()
        #@require Cbc="9961bab8-2fa3-5c5a-9d89-47fab24efd76"        include(joinpath(@__DIR__, "subsolvers", "cbc.jl"))
        @require Clp="e2554f3b-3117-50c0-817c-e040a3ddf72d"        include(joinpath(@__DIR__, "subsolvers", "clp.jl"))
        @require CPLEX="a076750e-1247-5638-91d2-ce28b192dca0"      include(joinpath(@__DIR__, "subsolvers", "cplex.jl"))
        @require GLPK="60bf3e95-4087-53dc-ae20-288a0d20c6a6"       include(joinpath(@__DIR__, "subsolvers", "glpk.jl"))
        @require Gurobi="2e9cd046-0924-5485-92f1-d5272153d98b"     include(joinpath(@__DIR__, "subsolvers", "gurobi.jl"))
        @require Hypatia="b99e6be6-89ff-11e8-14f8-45c827f4f8f2"    include(joinpath(@__DIR__, "subsolvers", "hypatia.jl"))
        @require KNITRO="67920dd8-b58e-52a8-8622-53c4cffbe346"     include(joinpath(@__DIR__, "subsolvers", "knitro.jl"))
        @require MosekTools="1ec41992-ff65-5c91-ac43-2df89e9693a4" include(joinpath(@__DIR__, "subsolvers", "mosek.jl"))
        @require Xpress="9e70acf3-d6c9-5be6-b5bd-4e2c73e3e054"     include(joinpath(@__DIR__, "subsolvers", "xpress.jl"))
    end

    #include("precompile.jl")
    #_precompile_()
end
