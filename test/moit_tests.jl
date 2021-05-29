
const unit_excludes = String[
    "number_threads",                 # EAGO won't support number of threads in near future
    "raw_status_string",              # TODO: ADD internal status states to EAGO
    "variablenames",                  # SEEMS LIKE A FALSE PASS

    "solve_qcp_edge_cases",           # Need to update nonconvex...
    "solve_qp_zero_offdiag",
     "solve_qp_edge_cases",

    "solve_zero_one_with_bounds_3",   # GLPK has a non-standard return code
    "solve_result_index",             # TODO: Should throw error when querying for multiple results... (expected behavior?)
]

const contlinear_excludes = String[
    "partial_start",                 # EAGO doesn't support VariablePrimalStart
]

const intlinear_excludes = String[
    "indicator1",    # doesn't currently support indicator sets
    "indicator2",    # can't check using Cbc until https://github.com/jump-dev/Cbc.jl/issues/151 is resolved
    "indicator3",
    "indicator4",
    "int2",
]

const contconic_excludes = String[
    "dualexp",  # Not directly bridged to common cones
    "dualpow",
    "logdet",
    "rootdet",
]

const contquadratic_excludes = String[
    "ncqcp",
    #"qp",
    "socp",
    "qcp",
]

const nlp_excludes = String[
    "feasibility_sense_with_objective_and_hessian",
    "feasibility_sense_with_no_objective_and_no_hessian",
    "feasibility_sense_with_objective_and_no_hessian",
    "hs071_no_hessian",
    "feasibility_sense_with_no_objective_and_with_hessian",
    "nlp_objective_and_moi_objective",
    "hs071",
    "hs071_hessian_vector_product_test"
]

function test_moi(T::Type{<:Real}; solver_options...)

    optimizer = MOIU.CachingOptimizer(MOIU.UniversalFallback(MOIU.Model{T}()), EAGO.Optimizer())
    MOI.set(optimizer, MOI.RawParameter(:verbosity), 0)

    tol = 2sqrt(sqrt(eps(T)))
    config = MOIT.TestConfig{T}(
    atol = tol,
    rtol = tol,
    solve = true,
    query = true,
    modify_lhs = false,
    duals = false,
    infeas_certificates = false,
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

    @testset "nlp" begin
        MOIT.nlptest(MOIB.full_bridge_optimizer(optimizer, T), config, nlp_excludes)
    end
end

# Test with mip_solver = Cbc as it supports SOS1 & SOS2 constraints
# TODO: Use bridges for SOS1 & SOS2 constraint if unsupported
# Need to test with GLPK as well to ensure subsolver supports constraint
# coefficient modification.
test_moi(Float64)
