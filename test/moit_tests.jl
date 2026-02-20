# Copyright (c) 2018 Matthew Wilhelm, Robert Gottlieb, Dimitri Alston,
# Alireza Miraliakbar, Matthew Stuber, and the University of Connecticut (UConn)
# This code is licensed under the MIT license (see LICENSE.md for full details).
################################################################################
# EAGO
# A development environment for robust and global optimization
# https://github.com/PSORLab/EAGO.jl
################################################################################
# test/moit_tests.jl
# The file runs all of the MathOptInterface tests on EAGO.
################################################################################

module TestEAGO

using Test

import EAGO
import MathOptInterface as MOI

"""
    runtests()

This function runs all functions in this Module starting with `test_`.
"""
function runtests()
    for name in names(@__MODULE__; all = true)
        if startswith("$(name)", "test_")
            @testset "$(name)" begin
                getfield(@__MODULE__, name)()
            end
        end
    end
    return
end

"""
    test_runtests()

This function runs all the tests in MathOptInterface.Test.

Pass arguments to `exclude` to skip tests for functionality that is not
implemented or that your solver doesn't support.
"""
function test_runtests()
    model = MOI.instantiate(
        EAGO.Optimizer;
        with_bridge_type = Float64,
        with_cache_type = Float64,
    )
    MOI.set(model, MOI.Silent(), true)
    MOI.Test.runtests(
        model,
        MOI.Test.Config(;
            atol = 1e-3,
            rtol = 1e-3,
            exclude = Any[
                MOI.DualObjectiveValue,
                MOI.ConstraintBasisStatus,
                MOI.VariableBasisStatus,
                MOI.ConstraintDual,
            ],
        );
        exclude = String[
            # Okay to exclude: returns INFEASIBLE_OR_UNBOUNDED instead of
            # INFEASIBLE
            "test_conic_NormInfinityCone_INFEASIBLE",
            "test_conic_NormOneCone_INFEASIBLE",
            # EAGO doesn't store result if the termination status is not FEASIBLE_POINT
            "test_linear_DUAL_INFEASIBLE",
            "test_linear_DUAL_INFEASIBLE_2",
            # Okay to exclude: returns INVALID_MODEL instead of INFEASIBLE
            "test_constraint_ZeroOne_bounds_3",
            # Okay to exclude: these tests throw the following warning and don't
            # terminate quickly.
            # ┌ Warning: At least one branching variable is unbounded. This will interfere with EAGO's global
            # │ optimization routine and may cause unexpected results. Bounds have been automatically
            # │ generated at +/- 1E10 for all unbounded variables, but tighter user-defined bounds are
            # │ highly recommended. To disable this warning and the automatic generation of bounds, use
            # │ the option `unbounded_check = false`.
            # └ @ EAGO ~/.julia/dev/EAGO/src/eago_optimizer/optimize/nonconvex/stack_management.jl:256
            "test_objective_qp_ObjectiveFunction_edge_cases",
            "test_quadratic_SecondOrderCone_basic",
            # TODO: wrong solutions. Likely a bug in EAGO.jl
            "test_conic_NormInfinityCone_VectorAffineFunction",
            "test_conic_NormInfinityCone_VectorOfVariables",
            "test_conic_NormOneCone_VectorAffineFunction",
            "test_conic_NormOneCone_VectorOfVariables",
            "test_conic_linear_VectorAffineFunction",
            "test_linear_Semicontinuous_integration",
            "test_linear_Semiinteger_integration",
            "test_linear_integer_solve_twice",
            "test_linear_integration",
            "test_quadratic_duplicate_terms",
            "test_quadratic_integration",
            "test_quadratic_nonhomogeneous",
            "test_modification_affine_deletion_edge_cases",
            # These test check `isapprox` between the `VectorNonlinearFunction`
            # and another one that was modified by bridges so they fail
            # See https://github.com/jump-dev/MathOptInterface.jl/issues/2553
            "test_basic_VectorNonlinearFunction_HyperRectangle",
            "test_basic_VectorNonlinearFunction_NormInfinityCone",
            "test_basic_VectorNonlinearFunction_NormOneCone",
            # TODO: bug related to unbounded_check
            "test_nonlinear_expression_quartic",
            # TODO: bug related to reform_epigraph_min
            "test_nonlinear_expression_overrides_objective",
            # Remove after MOI v1.31.3
            "test_nonlinear_expression_hs110",
            # TODO: wrong error thrown. Likely (trivial) bug in MOI wrapper
            "test_model_LowerBoundAlreadySet",
            "test_model_UpperBoundAlreadySet",
        ],
        exclude_tests_after = v"1.31.2",
    )
    return
end

end

@testset "MOI" begin
    TestEAGO.runtests()
end
