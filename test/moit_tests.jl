
const unit_excludes = String[
    "number_threads",                 # EAGO won't support number of threads in near future
    "raw_status_string",              # TODO: ADD internal status states to EAGO
    "solve_unbounded_model",          # CBC returns infeasible or unbounded for linear...

    "solve_zero_one_with_bounds_3",   # GLPK has a non-standard return code
    "solve_result_index",             # TODO: Should throw error when querying for multiple results... (expected behavior?)

    "solve_qcp_edge_cases",           # Not box constrained NLP type problems...
    "solve_qp_zero_offdiag",
    "solve_qp_edge_cases",

    # Passing
    #=
    "solve_objbound_edge_cases",
    "solve_integer_edge_cases",
    "solve_zero_one_with_bounds_1",
    "solve_zero_one_with_bounds_2",
    "solve_duplicate_terms_vector_affine",
    "solve_with_upperbound",
    "solve_single_variable_dual_max",
    "solve_single_variable_dual_min",
    "solve_farkas_lessthan",
    "solve_farkas_interval_upper",
    "solve_farkas_greaterthan",
    "solve_with_lowerbound",
    "solve_twice",
    "solve_farkas_variable_lessthan_max",
    "solve_farkas_variable_lessthan",
    "delete_soc_variables",
    "solve_farkas_equalto_lower",
    "solve_farkas_equalto_upper",
    "delete_nonnegative_variables",
    "solve_affine_greaterthan",
    "solve_farkas_interval_lower",
    "add_variable",
    "solve_singlevariable_obj",
    "solve_constant_obj",
    "solver_name",
    "delete_variable",
    "solve_time",
    "solve_duplicate_terms_obj",
    "delete_variables",
    "add_variables",
    "feasibility_sense",
    "solve_affine_equalto",
    "max_sense",
    "variablenames",
    "silent",
    "update_dimension_nonnegative_variables",
    "solve_blank_obj",
    "min_sense",
    "solve_affine_interval",
    "solve_duplicate_terms_scalar_affine",
    "time_limit_sec",
    "getvariable",
    "getconstraint",
    "get_objective_function",
    "solve_affine_lessthan",
    =#
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
