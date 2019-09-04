module EAGO

    import MathOptInterface

    using JuMP
    import JuMP._Derivatives: operators
    using Ipopt, GLPK

    using SparseArrays: SparseMatrixCSC, spzeros, rowvals, nzrange, nonzeros
    using StaticArrays: SVector
    using LinearAlgebra: eigmin

    import Calculus.symbolic_derivatives_1arg
    import Printf.@sprintf
    import Reexport.@reexport

    import IntervalArithmetic: +, -, *, /, convert, in, isempty, one, zero,
                               real, eps, max, min, abs, exp,
                               expm1, log, log2, log10, log1p, sqrt,
                               sin, cos, tan, min, max, sec, csc, cot, step,
                               sign, dist, mid, pow, Interval, sinh, cosh, ∩,
                               IntervalBox, bisect, isdisjoint, ^, exp2, exp10,
                               tanh, asinh, cosh, atanh


    const MOI = MathOptInterface
    const MOIU = MOI.Utilities

    include("mccormick_library/mccormick.jl")
    using .McCormick

    import Base: eltype, copy, length

    # Add storage types for EAGO optimizers
    export NodeBB, get_history, get_lower_bound, get_upper_bound, get_lower_time,
           get_upper_time, get_preprocess_time, get_postprocess_time, get_lower_bound, get_solution_time,
           get_iteration_number, get_node_count, get_absolute_gap, get_relative_gap
    include("eago_optimizer/branch_bound/node_bb.jl")
    include("eago_optimizer/branch_bound/subproblem_info.jl")
    include("eago_relaxations/relax_scheme.jl")
    include("eago_optimizer/moi_wrapper/optimizer.jl")

    # Routines for (branch-and-bound/cut)
    include("eago_optimizer/branch_bound/branch_bound.jl")

    # Routines for loading relax models
    include("eago_relaxations/relax_model.jl")

    # MOI wrappers and basic optimization routines
    include("eago_optimizer/moi_wrapper/moi_wrapper.jl")

    # Domain reduction subroutines
    include("eago_optimizer/domain_reduction/domain_reduction.jl")         # special character warning

    # Default subroutines for optimizers
    include("eago_optimizer/default_optimizer/default_optimizer.jl")

    # Adds the parametric interval methods
    export parametric_interval_params, param_intv_contractor
    include("parametric_interval/parametric_interval.jl")

    # Solve for implicit function optimization
    export ImplicitLowerEvaluator, build_lower_evaluator!,
           ImplicitUpperEvaluator, MidPointUpperEvaluator,
           build_implicitupperevaluator!, build_midpointupperevaluator!, solve_implicit
    include("eago_optimizer/implicit_optimizer/implicit_optimizer.jl")

    # Import the script solving utilities
    include("eago_script/script.jl")


    # Routines for solving SIPs
    export SIP_Options, SIP_Result, explicit_sip_solve, implicit_sip_solve
    include("eago_semiinfinite/semi_infinite.jl")
end
