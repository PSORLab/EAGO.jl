#=
const unit_excludes = String[
    "number_threads",                 # EAGO won't support number of threads in near future
    "raw_status_string",              # TODO: ADD internal status states to EAGO
    "solve_unbounded_model",          # CBC returns infeasible or unbounded for linear...
    "solve_zero_one_with_bounds_3",   # GLPK has a non-standard return code
    "solve_result_index",             # TODO: Should throw error when querying for multiple results... (expected behavior?)
    "solve_qcp_edge_cases",           # Not box constrained NLP type problems...
    "solve_qp_zero_offdiag",
    "solve_qp_edge_cases",
    "solve_affine_deletion_edge_cases"   # TODO: Fix this
]

const contlinear_excludes = String[
    "partial_start",  # EAGO doesn't support VariablePrimalStart
    "linear1",   # TODO: Fix this
 
    #=
    "linear13",
    "linear8a",
    "linear14",
    "linear6",
    "linear4",
    "linear3",
    "linear9",
    "linear8c",
    "linear2",
    "linear12",
    "linear7",
    "linear8b",
    "linear10b",
    "linear10",
    "linear15",
    "linear5",
    "linear11"
    =#
]

const intlinear_excludes = String[
    "indicator1",  # doesn't currently support indicator sets
    "indicator2",  # can't check using Cbc until https://github.com/jump-dev/Cbc.jl/issues/151 is resolved
    "indicator3",
    "indicator4",

    "int2",        # currently doesn't support sos1 or sos2 constraints
    "int3",   # TODO: Fix this

    # Passing
    #"semiinttest",
    #"semiconttest",
    #"int1",
    #"knapsack"
]

const contconic_excludes = String[
    "dualexp",  # Not directly bridged to common cones
    "dualpow",
    "logdet",
    "rootdet",
    "sdp",
    "normnuc",
    "exp",
    "soc",
    "normspec",
    "relentr",
    "rsoc",
    "pow",
    "geomean"
]

const contquadratic_excludes = String[
    "ncqcp",
    "qp",
    "socp",
    "qcp",
]

function test_moi(T::Type{<:Real}; solver_options...)

    optimizer = MOIU.CachingOptimizer(MOIU.UniversalFallback(MOIU.Model{T}()), EAGO.Optimizer())
    MOI.set(optimizer, MOI.RawParameter("verbosity"), 0)

    tol = 2sqrt(sqrt(eps(T)))
    config = MOIT.Config(T;
    atol = tol,
    rtol = tol,
    #solve = true,
    #query = true,
    #modify_lhs = false,
    #duals = false,
    #infeas_certificates = false,
    )

    @testset "unit tests" begin
        MOIT.unittest(MOIB.full_bridge_optimizer(optimizer, T), config, unit_excludes)
    end
    
    @testset "continuous linear tests" begin
        MOIT.contlineartest(MOIB.full_bridge_optimizer(optimizer, T), config, contlinear_excludes)
    end

    @testset "mixed-integer linear tests" begin
        MOIT.intlineartest(MOIB.full_bridge_optimizer(optimizer, T), config, intlinear_excludes)
    end
 
    @testset "continuous conic tests" begin
        MOIT.contconictest(MOIB.full_bridge_optimizer(optimizer, T), config, contconic_excludes)
    end

    @testset "continuous quadratic tests" begin
        MOIT.contquadratictest(MOIB.full_bridge_optimizer(optimizer, T), config, contquadratic_excludes)
    end
end

# Test with mip_solver = Cbc as it supports SOS1 & SOS2 constraints
# TODO: Use bridges for SOS1 & SOS2 constraint if unsupported
# Need to test with GLPK as well to ensure subsolver supports constraint
# coefficient modification.
test_moi(Float64)
=#

module TestEAGO

import EAGO
using MathOptInterface
using Test

const MOI = MathOptInterface
const OPTIMIZER = MOI.instantiate(MOI.OptimizerWithAttributes(EAGO.Optimizer, MOI.Silent() => true))
const BRIDGED = MOI.instantiate(MOI.OptimizerWithAttributes(EAGO.Optimizer, MOI.Silent() => true), with_bridge_type = Float64)
const CONFIG = MOI.Test.Config(atol = 1e-3, rtol = 1e-3, optimal_status = MOI.OPTIMAL, 
                               exclude = Any[MOI.DualObjectiveValue, MOI.ConstraintBasisStatus, MOI.VariableName, MOI.ConstraintName, MOI.delete,
                                             MOI.ConstraintDual, MOI.ListOfModelAttributesSet, MOI.add_constrained_variables])

"""
    runtests()

This function runs all functions in the this Module starting with `test_`.
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
                      exclude = [# IPOPT Inherited test exclusions
                                "test_model_ScalarFunctionConstantNotZero",
                                "test_solve_TerminationStatus_DUAL_INFEASIBLE",
                                "test_linear_VectorAffineFunction_empty_row",
                                "test_solve_DualStatus_INFEASIBILITY_CERTIFICATE_",
                                "test_model_LowerBoundAlreadySet",
                                "test_model_UpperBoundAlreadySet",
                                "test_model_copy_to_UnsupportedAttribute",
                                "test_model_copy_to_UnsupportedConstraint",
                                "test_objective_set_via_modify",
                                "test_model_ModelFilter_ListOfConstraintIndices",
                                "test_model_ModelFilter_ListOfConstraintTypesPresent",
                                # Cbc default test exlucisons
                                "_Indicator_",
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
                                "conic_NormOneCone_VectorAffineFunction",
                                "conic_NormOneCone_VectorOfVariables",
                                "conic_NormInfinityCone_VectorOfVariables",
                                "conic_NormInfinityCone_VectorAffineFunction",
                                "test_conic_NormInfinityCone_3",
                                "test_conic_NormOneCone",
                                "test_conic_linear_VectorOfVariables",
                                "test_conic_linear_VectorOfVariables_2",

                                "linear_integer_solve_twice",
                                "linear_integration",
                                "conic_linear_VectorAffineFunction",
                                "conic_linear_VectorAffineFunction_2",

                                "test_linear_integer_knapsack",
                                "test_linear_integer_integration",
                                "test_linear_inactive_bounds",
                                "test_linear_add_constraints",
                                "test_linear_FEASIBILITY_SENSE",

                                "test_constraint_ScalarAffineFunction_EqualTo",
                                "test_constraint_ScalarAffineFunction_GreaterThan",
                                "test_constraint_ScalarAffineFunction_Interval",
                                "test_constraint_ScalarAffineFunction_LessThan",
                                "test_constraint_ScalarAffineFunction_duplicate",
                                "test_constraint_VectorAffineFunction_duplicate",

                                "test_quadratic_SecondOrderCone_basic",
                                "test_quadratic_constraint_GreaterThan",
                                "test_quadratic_constraint_LessThan",
                                "test_quadratic_constraint_basic",
                                "test_quadratic_constraint_minimize",
                                "test_quadratic_duplicate_terms",
                                "test_quadratic_integration",
                                "test_quadratic_nonconvex_constraint_integration",
                                "test_quadratic_nonhomogeneous",
                                "test_quadratic_constraint_integration",

                                "test_variable_solve_Integer_with_lower_bound",
                                "test_variable_solve_Integer_with_upper_bound",
                                "test_variable_solve_with_lowerbound",
                                "test_variable_solve_with_upperbound",

                                "test_solve_result_index",

                                "test_objective_ObjectiveFunction_duplicate_terms",
                                "test_objective_ObjectiveFunction_constant",
                                "test_objective_ObjectiveFunction_VariableIndex",

                                "test_solve_ObjectiveBound_MAX_SENSE_IP",
                                "test_solve_ObjectiveBound_MAX_SENSE_LP",
                                "test_solve_ObjectiveBound_MIN_SENSE_IP",
                                "test_solve_ObjectiveBound_MIN_SENSE_LP",                                                                                          

                                "test_constraint_qcp_duplicate_diagonal",
                                "test_constraint_qcp_duplicate_off_diagonal",
                                "test_objective_get_ObjectiveFunction_ScalarAffineFunction",
                                "test_objective_qp_ObjectiveFunction_edge_cases",
                                "test_objective_qp_ObjectiveFunction_zero_ofdiag",
                                ],
                      exclude_tests_after = v"0.10.5")
end

end

@testset "MOI" begin
    TestEAGO.runtests()
end            