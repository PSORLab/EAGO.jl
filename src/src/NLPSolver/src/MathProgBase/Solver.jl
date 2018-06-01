"""
    EAGO_NLPSolver

Main solver type for EAGO global optimization. Contains all options that are not
modified over the course of the optimization problem. The fields are given below:
* `BnBSolver::BnBSolver`: The BnB solver object that that is modified then passed
                          to the solve function EAGOBranchBound. (Default = BnBSolver())
* `Implicit_Options::ImplicitSolver`: Solver options for implicit bounding routines. (Default = ImplicitSolver())
* `LBD_func_relax::String`: Relaxation type used in lower bounding problem. (Default = "NS-STD-OFF")
* `LBDsolvertype::String`: Type of problem relaxation to use when solving lower problem. (Default = "LP")
* `UBDsolvertype::String`: Type of problem relaxation to use when solving upper problem. (Default = "MPBNonlinear")
* `LP_solver`: LP solver for use in contraction routines. (Default = ClpSolver())
* `abs_tol_LBD::Float64`: Absolute tolerance spec for lower subproblem. (Default = 1E-5)
* `max_int_LBD::Int64`: Maximum iterations for lower subproblem. (Default = 5E5)
* `UBD_full_depth`: Depth below which problems are solved to feasilibity only (Default = 100)
* `abs_tol_UBD::Float64`: Absolute tolerance spec for upper subproblem. (Default = 1E-5)
* `max_int_UBD::Int64`: Maximum iterations for upper subproblem. (Default = 5E5)
* `STD_RR_depth::Int64`: Depth in tree to perform standard range reduction until. (Default = 1E10)
* `probe_depth::Int64`: Depth in tree to perform LP probing until. (Default = 3)
* `variable_depth::Int64`: Depth in tree to OBBT until. (Default = 1E15)
* `dual_tol::Float64`: Tolerance for recognizing a dual as on the bound. (Default = 1E-7)
* `DAG_depth::Int64`: Depth in tree to run DAG constraint propagation. (Default = 1E3)
* `DAG_pass::Int64`: Number of passes to run DAG constraint propagation. (Default = 3)
* `max_reduce_rept::Int64`: Maximum number of times to repeat tightening. (Not used currently.)
* `tol_reduce_rept::Float64`: Tolerance for repeating a node. (Not used currently.)
* `atol::Float64`: Absolute tolerance for termination. (Default = 1E-4)
* `rtol::Float64`: Relative tolerance for termination. (Default = 1E-4)
* `verbosity::String`: Verbosity of solution routine passed to BnB solve. (Default = "Normal")
* `iter_limit::Int64`: Iteration limit for branch and bound. (Default = "Normal")
* `node_limit::Int64`: Node limit for branch and bound. (Default = "Normal")
* `UBDsolver`: Default upper bounding solver
* `validated::Bool`: Flag indicating the interval calculation should be correctly rounded.
"""
type EAGO_NLPSolver <: AbstractMathProgSolver

    # Branch and Bound Solver Object
    BnBSolver::BnBSolver                  # The BnB solver object that that is modified then passed to the solve function EAGOBranchBound. (Default = BnBSolver())

    # Solver types
    LBD_func_relax::String                # Relaxation type used in lower bounding problem. (Default = "NS-STD-OFF")
    LBDsolvertype::String
    UBDsolvertype::String

    # lower bounding problem options
    LP_solver
    abs_tol_LBD::Float64
    max_int_LBD::Int64

    # upper bounding problem options
    UBD_full_depth::Int64
    abs_tol_UBD::Float64
    max_int_UBD::Int64

    # OBBT options
    STD_RR_depth::Int64
    probe_depth::Int64
    variable_depth::Int64
    dual_tol::Float64

    # FBBT options
    DAG_depth::Int64
    DAG_pass::Int64
    max_reduce_rept::Int64
    tol_reduce_rept::Float64

    # Limits & Tolerances
    atol::Float64
    rtol::Float64
    verbosity::String
    iter_limit::Int64
    node_limit::Int64

    UBDsolver
    validated::Bool

    ImplicitFlag::Bool
    PSmcOpt::mc_opts
    PIntOpt::PIntvParams

    Imp_RR_depth::Int64
    Imp_probe_depth::Int64

    SubGradRefine::Bool
end

function EAGO_NLPSolver(;

    BnBobject = BnBSolver(),

    LBD_func_relax = "NS-STD",
    LBDsolvertype = "LP",
    UBDsolvertype = "MPBNonlinear",

    LP_solver = ClpSolver(),
    abs_tol_LBD = 1E-5,
    max_int_LBD = 5E5,

    UBD_full_depth = 100,
    abs_tol_UBD = 1E-5,
    max_int_UBD = 5E5,

    STD_RR_depth = -1E10,
    probe_depth = -3,
    variable_depth = 1E15,
    dual_tol = 1E-7,

    DAG_depth = -10,
    DAG_pass = 3,
    max_reduce_rept = 5,
    tol_reduce_rept = 0.5,

    atol = 1E-4,
    rtol = 1E-4,
    verbosity = "Normal",
    iter_limit = 1E9,
    node_limit = 1E7,
    UBDsolver = IpoptSolver(print_level=0),
    validated = false,

    ImplicitFlag = false,
    PSmcOpt = mc_opts(),
    PIntOpt = PIntvParams(0,0),
    Imp_RR_depth = -1,
    Imp_probe_depth = -1,
    SubGradRefine = false)

    set_to_default!(BnBobject)
    set_Verbosity!(BnBobject,verbosity)

    BnBobject.BnB_atol = atol
    BnBobject.BnB_rtol = rtol
    BnBobject.iter_lim = true
    BnBobject.max_iter = iter_limit
    BnBobject.max_nodes = node_limit


    nlp_solver = EAGO_NLPSolver(deepcopy(BnBobject),
                                LBD_func_relax, LBDsolvertype, UBDsolvertype,
                                LP_solver, abs_tol_LBD, max_int_LBD, UBD_full_depth,
                                abs_tol_UBD, max_int_UBD, STD_RR_depth,
                                probe_depth, variable_depth, dual_tol,
                                DAG_depth, DAG_pass,
                                max_reduce_rept, tol_reduce_rept,
                                atol, rtol, verbosity,
                                iter_limit, node_limit, UBDsolver, validated,
                                ImplicitFlag, PSmcOpt, PIntOpt, Imp_RR_depth,
                                Imp_probe_depth, SubGradRefine)

    #Solver_Relax_Valid_LBD!(nlp_solver)

    return nlp_solver
end
