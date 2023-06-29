module TestEAGO

import EAGO
using MathOptInterface
using Test

const MOI = MathOptInterface
const OPTIMIZER = MOI.instantiate(MOI.OptimizerWithAttributes(EAGO.Optimizer, MOI.Silent() => true))
const BRIDGED = MOI.instantiate(MOI.OptimizerWithAttributes(EAGO.Optimizer, MOI.Silent() => true), with_bridge_type = Float64)
const CONFIG = MOI.Test.Config(atol = 1E-3, rtol = 1E-3, optimal_status = MOI.OPTIMAL, 
                               exclude = Any[MOI.DualObjectiveValue, MOI.ConstraintBasisStatus, MOI.ConstraintName, MOI.delete,
                                             MOI.ConstraintDual, MOI.ListOfModelAttributesSet, MOI.add_constrained_variables])

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
end

"""
    test_runtests()

This function runs all the tests in MathOptInterface.Test.

Pass arguments to `exclude` to skip tests for functionality that is not
implemented or that your solver doesn't support.
"""
function test_runtests()
    MOI.Test.runtests(BRIDGED, CONFIG, 
                      exclude = [# Ipopt inherited test exclusions
                                "test_model_ScalarFunctionConstantNotZero",
                                "test_solve_TerminationStatus_DUAL_INFEASIBLE",
                                "test_linear_VectorAffineFunction_empty_row",
                                "test_solve_DualStatus_INFEASIBILITY_CERTIFICATE_",
                                "test_model_LowerBoundAlreadySet",
                                "test_model_UpperBoundAlreadySet",
                                "test_model_copy_to_UnsupportedAttribute",
                                "test_model_copy_to_UnsupportedConstraint",
                                "test_objective_set_via_modify",                                
                                "test_conic_linear_VectorOfVariables_2",

                                # Cbc default test exclusions
                                "test_linear_Indicator_",
                                "test_linear_SOS1_integration",
                                "test_linear_SOS2_integration",
                                "test_solve_SOS2_add_and_delete",
                                "test_conic_NormInfinityCone_INFEASIBLE",
                                "test_conic_NormOneCone_INFEASIBLE",
                                "test_solve_TerminationStatus_DUAL_INFEASIBLE",

                                # EAGO test exclusions
                                "test_attribute_NumberOfThreads",
                                "test_modification_",
                                "test_linear_integration_delete_variables",

                                # EAGO exclusions to resolve (by adding conic support later and fixing twice solve issues)
                                "test_conic_NormOneCone_VectorAffineFunction",
                                "test_conic_NormOneCone_VectorOfVariables",
                                "test_conic_NormInfinityCone_VectorOfVariables",
                                "test_conic_NormInfinityCone_VectorAffineFunction",
                                "test_conic_NormInfinityCone_3",
                               r"^test_conic_NormOneCone$",
                               r"^test_conic_linear_VectorOfVariables$",
                               r"^test_conic_linear_VectorAffineFunction$",
                                "test_conic_linear_VectorAffineFunction_2",

                                "test_linear_integer_solve_twice",
                               r"^test_linear_integration$",

                                "test_quadratic_SecondOrderCone_basic",
                                "test_quadratic_constraint_GreaterThan",
                                "test_quadratic_constraint_LessThan",
                                "test_quadratic_constraint_minimize",
                                "test_quadratic_duplicate_terms",
                                "test_quadratic_integration",
                                "test_quadratic_nonconvex_constraint_integration",
                                "test_quadratic_homogeneous",
                                "test_quadratic_nonhomogeneous",
                                "test_quadratic_constraint_integration",
                                "test_quadratic_constraint_basic",
                                "test_quadratic_nonconvex_constraint_basic",                                                                               

                                "test_constraint_qcp_duplicate_diagonal",
                                "test_constraint_qcp_duplicate_off_diagonal",
                                "test_objective_get_ObjectiveFunction_ScalarAffineFunction",
                                "test_objective_qp_ObjectiveFunction_edge_cases",
                                "test_objective_qp_ObjectiveFunction_zero_ofdiag",

                                "test_model_ModelFilter_ListOfConstraintIndices",
                                "test_model_ModelFilter_ListOfConstraintTypesPresent",

                                # MOI constraint type exclusions
                                "test_cpsat_Circuit",
                                "test_cpsat_CountAtLeast",
                                "test_cpsat_Table",
                                "test_linear_Semicontinuous_integration",
                                "test_linear_Semiinteger_integration",

                                # EAGO B&B handles interval bounds internally
                                "test_linear_open_intervals",
                                "test_linear_variable_open_intervals",

                                # EAGO does not have constraint names
                                "test_model_Name_VariableName_ConstraintName"

                                ],
                      exclude_tests_after = v"1.18.0")
end

end

@testset "MOI" begin
    TestEAGO.runtests()
end            