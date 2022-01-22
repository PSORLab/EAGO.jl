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
                               exclude = Any[MOI.DualObjectiveValue, MOI.VariableName, MOI.DualObjectiveValue, MOI.delete,
                                             MOI.ConstraintFunction, MOI.ConstraintDual, MOI.ConstraintSet, 
                                             MOI.ListOfModelAttributesSet, MOI.ListOfConstraintIndices, MOI.ListOfConstraintTypesPresent,
                                             MOI.add_constrained_variables])

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
                                #  - Excluded because this test is optional
                                "test_model_ScalarFunctionConstantNotZero",
                                #  - Excluded because Ipopt returns NORM_LIMIT instead of
                                #    DUAL_INFEASIBLE
                                "test_solve_TerminationStatus_DUAL_INFEASIBLE",
                                #  - Excluded because Ipopt returns INVALID_MODEL instead of
                                #    LOCALLY_SOLVED
                                "test_linear_VectorAffineFunction_empty_row",
                                "test_solve_DualStatus_INFEASIBILITY_CERTIFICATE_",
                                #  - Excluded due to upstream issue
                                "test_model_LowerBoundAlreadySet",
                                "test_model_UpperBoundAlreadySet",
                                #  - CachingOptimizer does not throw if optimizer not attached
                                "test_model_copy_to_UnsupportedAttribute",
                                "test_model_copy_to_UnsupportedConstraint",
                                # EAGO test exclusions
                                "test_attribute_NumberOfThreads",
                                "test_modification_",
                                "test_quadratic_",
                                "test_variable_delete_",
                                "test_model_ModelFilter_ListOfConstraintIndices",
                                "test_model_ModelFilter_ListOfConstraintTypesPresent",
                                "test_add_constrained_variables_vector",
                                "test_basic_ScalarAffineFunction_EqualTo",
                                "test_basic_ScalarAffineFunction_GreaterThan",
                                "test_basic_ScalarAffineFunction_LessThan",
                                "test_basic_ScalarQuadraticFunction_EqualTo",
                                "test_basic_ScalarQuadraticFunction_GreaterThan",
                                "test_basic_ScalarQuadraticFunction_LessThan",
                                "test_basic_VariableIndex_EqualTo",
                                "test_basic_VariableIndex_GreaterThan",
                                "test_basic_VariableIndex_LessThan",
                                "test_conic_NormInfinityCone_3",
                                "test_conic_NormInfinityCone_VectorAffineFunction",
                                "test_conic_NormInfinityCone_VectorOfVariables",
                                "test_conic_NormOneCone",
                                "test_conic_NormOneCone_VectorAffineFunction",
                                "test_conic_NormOneCone_VectorOfVariables",
                                "test_conic_linear_VectorAffineFunction",
                                "test_conic_linear_VectorAffineFunction_2",
                                "test_conic_linear_VectorOfVariables",
                                "test_constraint_ZeroOne_bounds_3",
                                "test_constraint_qcp_duplicate_diagonal",
                                "test_constraint_qcp_duplicate_off_diagonal",
                                "test_linear_Interval_inactive",
                                "test_linear_Semicontinuous_integration",
                                "test_linear_Semiinteger_integration",
                                "test_linear_add_constraints",
                                "test_linear_inactive_bounds",
                                "test_linear_integer_integration",
                                "test_linear_integer_knapsack",
                                "test_linear_integer_solve_twice",
                                "test_linear_integration",
                                "test_linear_integration_2",
                                "test_linear_integration_Interval",
                                "test_linear_integration_modification",
                                "test_linear_transform",
                                "test_model_empty",
                                "test_model_is_valid",
                                "test_objective_get_ObjectiveFunction_ScalarAffineFunction",
                                "test_objective_qp_ObjectiveFunction_edge_cases",
                                "test_objective_qp_ObjectiveFunction_zero_ofdiag",
                                "test_objective_set_via_modify",
                                ],
                      exclude_tests_after = v"0.10.5")
end

end

@testset "MOI" begin
    TestEAGO.runtests()
end            