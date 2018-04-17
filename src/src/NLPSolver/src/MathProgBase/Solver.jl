"""
    ImplicitSolver

ImplicitSolver is an option storage type utilized by the EAGO_NLPSolver to set
the options used for solve NLPs with embedded implicit calculations. By default,
it's flag field is set to false and the other fields are effectively unset. It's fields
are as follows:
* `flag::Bool`: A flag indicating whether or not implicit bounding routines should
                be used to solve subproblems.
* `opts::mc_opts`: Option type used for McCormick fixed point bounding routine.
* `ParamInt::Vector`: Options used for parametric interval calculations
                      (iteration #, equality tolerance, extended Newton parameter)
* `f`: The objective function used in implicit calculations. Takes the form: f(x,p)
       where x = y[1:nx] and p = y[(nx+1):end].
* `g`: Inequality constraints used in implicit calculations. Takes the form: g(x,p)
       where x = y[1:nx] and p = y[(nx+1):end]. Should be in the same order as
       constraints the explicit calculations omitting the constraints used to
       form the implicit bounds.
* `h`: Equality constraints used in implicit calculations to form relaxations of
       an implicit function. Takes the form: h(x,p) where x = y[1:nx] and
       p = y[(nx+1):end]. Should be in the same order as constraints the explicit
       calculations.
* `hj`: Jacobian of h w.r.t the state space variables 'z' and takes the form
* `nx`: The number of variables in the state space.
* `Intv_Cntr`: Style of interval contractor used in pre-processing. Options are
               "NewtonGS" for Newton Gauss-Siedel and "KrawczykCW" for a componentwise
               Krawczyk calculation.
* `Imp_RR_depth`: The depth in the tree to which implicit range reduction in which
                  the equality constraints are relaxed and g(x(p),p)<0, f(x(p),p)
                  forms the problem for range-reduction and only p bounds are tightened.
* `Imp_probe_depth`: The depth in the tree to which implicit range reduction in which
                      the equality constraints are relaxed and g(x(p),p)<0, f(x(p),p)
                      forms the problem for probing and only p bounds are tightened.
* `numConstr`: Number of constraints that remain after removing those which define
               the implicit function.
* `gL_Loc`: Indices at which a constraint has finite lower bounds.
* `gU_Loc`: Indices at which a constraint has finite upper bounds.
* `gL`: Lower bounds for inequality constraints.
* `gU`: Upper bounds for inequality constraints.
"""
type ImplicitSolver
    flag::Bool
    opts::mc_opts
    ParamInt::PIntvParams{Float64}
    f::Function
    g
    h
    hj
    nx
    Intv_Cntr
    Imp_RR_depth
    Imp_probe_depth
    numConstr
    gL_Loc
    gU_Loc
    gL
    gU
end
ImplicitSolver() = ImplicitSolver(false,mc_opts(),PIntvParams(0,0),
                                    x->x,x->x,x->x,x->x,0,
                                    "NewtonGS",-1,-1,[],[],[],[],[])

"""
    EAGO_NLPSolver

Main solver type for EAGO global optimization. Contains all options that are not
modified over the course of the optimization problem. The fields are given below:
* `BnBSolver::BnBSolver`: The BnB solver object that that is modified then passed
                          to the solve function EAGOBranchBound. (Default = BnBSolver())
* `Implicit_Options::ImplicitSolver`: Solver options for implicit bounding routines. (Default = ImplicitSolver())
* `LBD_func_relax::String`: Relaxation type used in lower bounding problem. (Default = "NS-STD-OFF")
* `LBD_problem_relax::String`: Type of problem relaxation to use when solving lower problem. (Default = "LP")
* `LBD_problem_solver::String`: Solver for use in problem lower problem. (Default = "Clp")
* `UBD_func_relax::String`: Relaxation type used in lower bounding problem. (Default = "Original")
* `UBD_problem_relax::String`: Type of problem relaxation to use when solving lower problem. (Default = "NLP2")
* `UBD_problem_solver::String`: Solver for use in problem upper problem. (Default = "Ipopt")
* `LP_solver`: LP solver for use in contraction routines. (Default = ClpSolver())
* `abs_tol_LBD::Float64`: Absolute tolerance spec for lower subproblem. (Default = 1E-5)
* `max_int_LBD::Int64`: Maximum iterations for lower subproblem. (Default = 5E5)
* `UBD_feas_depth`: Depth below which problems are solved to feasilibity only (Default = 100)
* `abs_tol_UBD::Float64`: Absolute tolerance spec for upper subproblem. (Default = 1E-5)
* `max_int_UBD::Int64`: Maximum iterations for upper subproblem. (Default = 5E5)
* `STD_RR_depth::Int64`: Depth in tree to perform standard range reduction until. (Default = 1E10)
* `probe_depth::Int64`: Depth in tree to perform LP probing until. (Default = 3)
* `variable_depth::Int64`: Depth in tree to OBBT until. (Default = 1E15)
* `dual_tol::Float64: Tolerance for recognizing a dual as on the bound. (Default = 1E-7)
* `DAG_depth::Int64`: Depth in tree to run DAG constraint propagation. (Default = 1E3)
* `DAG_pass::Int64`: Number of passes to run DAG constraint propagation. (Default = 3)
* `max_reduce_rept::Int64`: Maximum number of times to repeat tightening. (Not used currently.)
* `tol_reduce_rept::Float64`: Tolerance for repeating a node. (Not used currently.)
* `atol::Float64`: Absolute tolerance for termination. (Default = 1E-4)
* `rtol::Float64`: Relative tolerance for termination. (Default = 1E-4)
* `verbosity::String`: Verbosity of solution routine passed to BnB solve. (Default = "Normal")
* `iter_limit::Int64`: Iteration limit for branch and bound. (Default = "Normal")
* `node_limit::Int64`: Node limit for branch and bound. (Default = "Normal")
"""
type EAGO_NLPSolver <: AbstractMathProgSolver

    # Branch and Bound Solver Object
    BnBSolver::BnBSolver
    Implicit_Options::ImplicitSolver

    # Solver types
    LBD_func_relax::String
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
    DAG_depth::Int64 # Depth in tree to run DAG constraint propagation
    DAG_pass::Int64  # Number of passes to run DAG constraint propagation
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
end

function EAGO_NLPSolver(;

    BnBobject = BnBSolver(),
    ImplicitOpts = ImplicitSolver(),

    LBD_func_relax = "NS-STD-OFF",
    LBDsolvertype = "LP",
    UBDsolvertype = "MPBNonlinear",

    LP_solver = ClpSolver(),
    abs_tol_LBD = 1E-5,
    max_int_LBD = 5E5,

    UBD_full_depth = 100,
    abs_tol_UBD = 1E-5,
    max_int_UBD = 5E5,

    STD_RR_depth = 1E10,
    probe_depth = 3,
    variable_depth = 1E15,
    dual_tol = 1E-7,

    DAG_depth = 1E3,
    DAG_pass = 3,
    max_reduce_rept = 5,
    tol_reduce_rept = 0.5,

    atol = 1E-4,
    rtol = 1E-4,
    verbosity = "Normal",
    iter_limit = 1E9,
    node_limit = 1E7,
    UBDsolver = IpoptSolver(print_level=0),
    validated = false)

    set_to_default!(BnBobject)
    set_Verbosity!(BnBobject,verbosity)

    BnBobject.BnB_atol = atol
    BnBobject.BnB_rtol = rtol

    nlp_solver = EAGO_NLPSolver(deepcopy(BnBobject),deepcopy(ImplicitOpts),
                                LBD_func_relax, LBDsolvertype, UBDsolvertype,
                                LP_solver, abs_tol_LBD, max_int_LBD, UBD_full_depth,
                                abs_tol_UBD, max_int_UBD, STD_RR_depth,
                                probe_depth, variable_depth, dual_tol,
                                DAG_depth, DAG_pass,
                                max_reduce_rept, tol_reduce_rept,
                                atol, rtol, verbosity,
                                iter_limit, node_limit, UBDsolver, validated)

    #Solver_Relax_Valid_LBD!(nlp_solver)
    #Solver_Relax_Valid_UBD!(nlp_solver)

    return nlp_solver
end
