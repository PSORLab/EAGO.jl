var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Introduction",
    "title": "Introduction",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#EAGO-Easy-Advanced-Global-Optimization-in-Julia-1",
    "page": "Introduction",
    "title": "EAGO - Easy Advanced Global Optimization in Julia",
    "category": "section",
    "text": ""
},

{
    "location": "index.html#Authors-1",
    "page": "Introduction",
    "title": "Authors",
    "category": "section",
    "text": "Matthew Wilhelm, Department of Chemical and Biomolecular Engineering,  University of Connecticut (UCONN)"
},

{
    "location": "index.html#Overview-1",
    "page": "Introduction",
    "title": "Overview",
    "category": "section",
    "text": "EAGO is a global and robust optimization platform based on McCormick relaxations. It contains the first widely accessible global optimization routine based on generalized McCormick relaxations. With the exception of calls to local solvers and linear algebra routines, EAGO is written entirely in native Julia. The solver is quite flexibly arranged so the end user can easily customize low-level routines."
},

{
    "location": "index.html#Installing-EAGO-1",
    "page": "Introduction",
    "title": "Installing EAGO",
    "category": "section",
    "text": "EAGO is registered Julia package and can be installed by running:julia> Pkg.add(\"EAGO\")"
},

{
    "location": "EAGOSolver/starting.html#",
    "page": "Getting Started",
    "title": "Getting Started",
    "category": "page",
    "text": ""
},

{
    "location": "EAGOSolver/starting.html#Solvers-with-a-guarantee-of-global-optimality-1",
    "page": "Getting Started",
    "title": "Solvers with a guarantee of global optimality",
    "category": "section",
    "text": "For unconstrained problems, the following upper bounding problem modes will provide a solution that is globally optimal to within an epsilon tolerance:Interval\nLP (Explicit)The following upper bounding problems solvers can also provide a valid solution to constrained problems:SNOPT\nIpopt\nMPBNonlinearThe following lower bounding solvers are currently under-construction and will likely furnish incorrect answers/errors.AlphaBB\nQuadratic\nIpoptLower bounding problem options are:Interval\nSNOPT\nLPThe solver contains a hook into JuMP that can be used to solve simple explicit problems as shown below: "
},

{
    "location": "EAGOSolver/starting.html#Solving-a-basic-problem-problem-1",
    "page": "Getting Started",
    "title": "Solving a basic problem problem",
    "category": "section",
    "text": "jumpmodel4 = Model(solver=EAGO_NLPSolver(LBD_func_relax = \"NS-STD-OFF\",\n                                         LBDsolvertype = \"LP\",\n                                         probe_depth = -1,\n                                         variable_depth = 1000,\n                                         DAG_depth = -1,\n                                         STD_RR_depth = -1))\n@variable(jumpmodel4, -200 <= x <= -100)\n@variable(jumpmodel4, 200 <= y <= 400)\n@constraint(jumpmodel4, -500 <= x+2y <= 400)\n@NLobjective(jumpmodel4, Min, x*y)\nstatus2 = solve(jumpmodel4)"
},

{
    "location": "EAGOSolver/SolverOpts.html#EAGO.EAGO_NLPSolver",
    "page": "Setting Solver Options",
    "title": "EAGO.EAGO_NLPSolver",
    "category": "type",
    "text": "EAGO_NLPSolver\n\nMain solver type for EAGO global optimization. Contains all options that are not modified over the course of the optimization problem. The fields are given below:\n\nBnBSolver::BnBSolver: The BnB solver object that that is modified then passed                         to the solve function EAGOBranchBound. (Default = BnBSolver())\nImplicit_Options::ImplicitSolver: Solver options for implicit bounding routines. (Default = ImplicitSolver())\nLBD_func_relax::String: Relaxation type used in lower bounding problem. (Default = \"NS-STD-OFF\")\nLBDsolvertype::String: Type of problem relaxation to use when solving lower problem. (Default = \"LP\")\nUBDsolvertype::String: Type of problem relaxation to use when solving upper problem. (Default = \"MPBNonlinear\")\nLP_solver: LP solver for use in contraction routines. (Default = ClpSolver())\nabs_tol_LBD::Float64: Absolute tolerance spec for lower subproblem. (Default = 1E-5)\nmax_int_LBD::Int64: Maximum iterations for lower subproblem. (Default = 5E5)\nUBD_full_depth: Depth below which problems are solved to feasilibity only (Default = 100)\nabs_tol_UBD::Float64: Absolute tolerance spec for upper subproblem. (Default = 1E-5)\nmax_int_UBD::Int64: Maximum iterations for upper subproblem. (Default = 5E5)\nSTD_RR_depth::Int64: Depth in tree to perform standard range reduction until. (Default = 1E10)\nprobe_depth::Int64: Depth in tree to perform LP probing until. (Default = 3)\nvariable_depth::Int64: Depth in tree to OBBT until. (Default = 1E15)\ndual_tol::Float64: Tolerance for recognizing a dual as on the bound. (Default = 1E-7)\nDAG_depth::Int64: Depth in tree to run DAG constraint propagation. (Default = 1E3)\nDAG_pass::Int64: Number of passes to run DAG constraint propagation. (Default = 3)\nmax_reduce_rept::Int64: Maximum number of times to repeat tightening. (Not used currently.)\ntol_reduce_rept::Float64: Tolerance for repeating a node. (Not used currently.)\natol::Float64: Absolute tolerance for termination. (Default = 1E-4)\nrtol::Float64: Relative tolerance for termination. (Default = 1E-4)\nverbosity::String: Verbosity of solution routine passed to BnB solve. (Default = \"Normal\")\niter_limit::Int64: Iteration limit for branch and bound. (Default = \"Normal\")\nnode_limit::Int64: Node limit for branch and bound. (Default = \"Normal\")\nUBDsolver: Default upper bounding solver\nvalidated::Bool: Flag indicating the interval calculation should be correctly rounded.\n\n\n\n"
},

{
    "location": "EAGOSolver/SolverOpts.html#",
    "page": "Setting Solver Options",
    "title": "Setting Solver Options",
    "category": "page",
    "text": "The majority of the basic solver settings are controlled via keyword inputs to the solver object. The EAGO_NLPSolver accepts the following options:EAGO_NLPSolver"
},

{
    "location": "EAGOSolver/MPB.html#",
    "page": "MathProgBase Inferace",
    "title": "MathProgBase Inferace",
    "category": "page",
    "text": ""
},

{
    "location": "EAGOSolver/MPB.html#Using-the-MathProgBase-Interface-1",
    "page": "MathProgBase Inferace",
    "title": "Using the MathProgBase Interface",
    "category": "section",
    "text": "The MathProgBase interface is recommended for optimizing nonsmooth objects and problems that can\'t be easily represented in JuMP. A EAGO model object should be first created from as solverf(x) = (x[1]-5)^2 + (x[2]-3)^2\ng(x) = [x[1] - x[2]]\nm1 = MathProgBase.NonlinearModel(EAGO_NLPSolver())\nMathProgBase.loadproblem!(m1, 2, 1, [0.0, 0.0], [10.0, 10.0],\n            [0.0], [1.0], :Min, f, g)\nMathProgBase.optimize!(m1)"
},

{
    "location": "EAGOSolver/MPB.html#EAGO.EAGO_NLP_Model",
    "page": "MathProgBase Inferace",
    "title": "EAGO.EAGO_NLP_Model",
    "category": "type",
    "text": "EAGO_NLP_Model\n\nThe model type to interface with JuMP and MathProgBase. This has fields:\n\nBnBModel::BnBModel: which stores information from the Branch and Bound Problem.\nOpts::EAGO_Inner_NLP: storage for various problem descriptors\nd: The storage for an AbstractNLPEvaluator\nstatus::Symbol: descriptor for model status (e.g. :optimal)\n\n\n\n"
},

{
    "location": "EAGOSolver/MPB.html#EAGO.EAGO_Inner_NLP",
    "page": "MathProgBase Inferace",
    "title": "EAGO.EAGO_Inner_NLP",
    "category": "type",
    "text": "EAGO_Inner_NLP\n\nThis is the inner storage object for the model interface to JuMP and MathProgBase. Effectively empty when initialized. It has the following fields:\n\ngL::Vector{Float64}: Lower constraint bounds\ngU::Vector{Float64}: Upper constraint bounds\ngexp::Int64: Expression for constraints\ngL_loc::Vector{Int64}: Index at which constraints have finite lower bounds\ngU_loc::Vector{Int64}: Index at which constraints have finite upper bounds\nnumVar::Int64: Number of variables in full problem.\nnumConstr::Int64: Number of constraints in full problem.\nvartypes::Vector{Symbol}: Type of variable array (currently not used)\nsense::Symbol: Minimization (:Min) or maximization (:Max)\nobj::Expr: Expression for objective provide by AbstractNLPEvaluator\nconstrs::Vector{Expr}: Expressions for constraints provide by AbstractNLPEvaluator\nDAG_tlist::TapeList: Tape list for DAG propagation on constraints\nf::Function: Objective function\ng::Function: Constraint function\nsolver::EAGO_NLPSolver: Solve storage object\n\n\n\n"
},

{
    "location": "EAGOSolver/MPB.html#EAGO-Model-Object-1",
    "page": "MathProgBase Inferace",
    "title": "EAGO Model Object",
    "category": "section",
    "text": "EAGO_NLP_ModelEAGO_Inner_NLP"
},

{
    "location": "EAGOSolver/MPB.html#MathProgBase.SolverInterface.optimize!-Tuple{EAGO.EAGO_NLP_Model}",
    "page": "MathProgBase Inferace",
    "title": "MathProgBase.SolverInterface.optimize!",
    "category": "method",
    "text": "MathProgBase.optimize!(s::EAGO_NLP_Model)\n\nOptimizes the s::EAGO_NLP_Model. May print console outputs depending on settings and other solution information is accessible via start solver interface functions.\n\n\n\n"
},

{
    "location": "EAGOSolver/MPB.html#Setup-API-for-MathProgBase-Interface-1",
    "page": "MathProgBase Inferace",
    "title": "Setup API for MathProgBase Interface",
    "category": "section",
    "text": "    MathProgBase.optimize!(s::EAGO_NLP_Model)"
},

{
    "location": "EAGOSolver/MPB.html#Access-and-Manipulation-functions-for-MathProgBase-Interface-1",
    "page": "MathProgBase Inferace",
    "title": "Access and Manipulation functions for MathProgBase Interface",
    "category": "section",
    "text": "Standard MathProgBase access functions (e.g. MathProgBase.getsolution) are extended by EAGO."
},

{
    "location": "EAGOSolver/SIP.html#",
    "page": "SIP Solver Interface",
    "title": "SIP Solver Interface",
    "category": "page",
    "text": ""
},

{
    "location": "EAGOSolver/SIP.html#EAGO.Explicit_SIP_Solve",
    "page": "SIP Solver Interface",
    "title": "EAGO.Explicit_SIP_Solve",
    "category": "function",
    "text": "Explicit_SIP_Solve\n\nSolves a semi-infinite program via the algorithm presented in Mitsos2011 using the EAGOGlobalSolver to solve the lower bounding problem, lower level problem, and the upper bounding problem. The options for the algorithm and the global solvers utilized are set by manipulating a SIPopt containing the options info. Inputs:\n\nf::Function: Objective in the decision variable. Takes a single argument                vector that must be untyped.\ngSIP::Function: The semi-infinite constraint. Takes two arguments: the first                   being a vector containing the decision variable and the                   second being a vector containing the uncertainity                   variables. The function must be untyped.\nX::Vector{Interval}: Box constraints for decision variables\nP::Vector{Interval}: Box constraints for uncertainty variables\nSIPopt::SIP_opts: Option type containing problem information\n\nReturns: A SIP_result composite type containing solution information.\n\n\n\n"
},

{
    "location": "EAGOSolver/SIP.html#Solving-Semi-Infinite-Programs-1",
    "page": "SIP Solver Interface",
    "title": "Solving Semi-Infinite Programs",
    "category": "section",
    "text": "Explicit_SIP_Solve"
},

{
    "location": "EAGOSolver/SIP.html#Example-of-solving-a-SIP-without-equality-constraints-1",
    "page": "SIP Solver Interface",
    "title": "Example of solving a SIP without equality constraints",
    "category": "section",
    "text": "using EAGO\n\nSIPopt1 = SIP_opts()\nsep1lu = EAGO_NLPSolver(probe_depth = -1,\n                        variable_depth = 1000,\n                        DAG_depth = -1,\n                        STD_RR_depth = -1,\n                        UBDsolvertype= \"Ipopt\")\nsep1lu.BnBSolver.Verbosity = \"Full\"\nsep1in = EAGO_NLPSolver(probe_depth = -1,\n                        variable_depth = 1000,\n                        DAG_depth = -1,\n                        STD_RR_depth = -1,\n                        UBDsolvertype= \"Ipopt\")\nsep1in.BnBSolver.Verbosity = \"Full\"\nSIPopt1.LLP_Opt = sep1in\nSIPopt1.LBP_Opt = sep1lu\nSIPopt1.UBP_Opt = sep1lu\nf1(x) = (1/3)*x[1]^2 + x[2]^2 + x[1]/2\ngSIP1(x,p) = (1.0-(x[1]^2)*(p[1]^2))^2 - x[1]*p[1]^2 - x[2]^2 + x[2]\nX1 = [MCInterval(-1000.0,1000.0),MCInterval(-1000.0,1000.0)]\nP1 = [MCInterval(0.0,1.0)]\nSIPoutput1 = Explicit_SIP_Solve(f1,gSIP1,X1,P1,SIPopt1)"
},

{
    "location": "EAGOSolver/SIP.html#Solving-Semi-Infinite-Programs-with-Equality-Constraints-1",
    "page": "SIP Solver Interface",
    "title": "Solving Semi-Infinite Programs with Equality Constraints",
    "category": "section",
    "text": "using EAGO\n\n# create the SIP option object for the solver\nSIPopt1 = SIP_opts()\n\n# create solver with specified options options for lower level problem\nsep1in = EAGO_NLPSolver(LBD_func_relax = \"NS-STD\",  # use standard McCormick relaxations\n                        LBDsolvertype = \"LP\",           # use an LP problem structure for relaxed problems\n                        UBDsolvertype = \"Ipopt\",        # use NLP solver upper bounds (currently preferred solver)\n                        probe_depth = -1,               # disable probing\n                        variable_depth = 1000,          # use duality based range reduction to a depth of 1000 (use to high depth recommended)\n                        DAG_depth = -1,                 # don\'t use a DAG contractor (I need to update this for implicit SIP)\n                        STD_RR_depth = -1,              # don\'t use standard range reduction (problems get quite large)\n                        verbosity = \"None\",             # specify printing level for global optimization problem\n                        validated = true,               # use numerically validated intervals\n                        atol = 1E-7,                    # absolute tolerance (May need to play with this)\n                        rtol = 1E-5)                    # relative tolerance (May need to play with this)\n\n# create a solver for the lower/upper problems\nsep1lu = EAGO_NLPSolver(LBD_func_relax = \"NS-STD\",\n                        LBDsolvertype = \"LP\",\n                        UBDsolvertype = \"Ipopt\",\n                        probe_depth = -1,\n                        variable_depth = 1000,\n                        DAG_depth = -1,\n                        STD_RR_depth = -1,\n                        verbosity = \"None\",\n                        validated = true,\n                        atol = 1E-7,\n                        rtol = 1E-5)\n\nSIPopt1.LLP_Opt = sep1in        # Set solver for use in lower level problem\nSIPopt1.LBP_Opt = sep1lu        # Set solver for use in lower bounding problem\nSIPopt1.UBP_Opt = sep1lu        # Set solver for use in upper bounding problem\n\nSIPopt1.eps_g0 = 0.9\nSIPopt1.tol = 1E-2              # SIP tolerance\nSIPopt1.r0 = 2.0                # reduction factor for SIP routine\nSIPopt1.kmax = 5                # maximum number of iteration for SIP routine\nSIPopt1.inn_tol = 0.0           # tolerance factor usually set to tolerance of inner program\n\n# 1D Example (7.4.1 from thesis)\n# solution f = -15.8077 @ y = 2.95275\nf(x) = (x[1]-3.5)^4 - 5*(x[1]-3.5)^3 - 2*(x[1]-3.5)^2 + 15*(x[1]-3.5)\nfunction h(x,y,p)\n    [y[1]-(x[1]-(x[1]^3)/6+(x[1]^5)/120)/sqrt(y[1])-p[1]]\nend\nfunction hj(x,y,p)\n    [1.0+(x[1]-(x[1]^3)/6+(x[1]^5)/120)/(2.0*sqrt(y[1]^3))]\nend\ngSIP(x,y,p) = y[1] + cos(x[1]-p[1]/90) - p[1]\nxBnds = [Interval(0.5,8.0)]\nyBnds = [Interval(68.8,149.9)]\npBnds = [Interval(80,120)]\nimpout1 = Implicit_SIP_Solve(f,h,hj,gSIP,xBnds,yBnds,pBnds,SIPopt1)\n\n# get solution values\nk = impout1.k                   # number of iterations\nUBD = impout1.UBD               # upper bound\nLBD = impout1.LBD               # lower bound\nfeas = impout1.feas             # is problem feasible?\nLBP_time = impout1.LBP_time     # time spent solving lower bounding problem\nLLP_time = impout1.LLP_time     # time spent solving lower level problem\nUBP_time = impout1.UBP_time     # time spent solving upper bounding problem\nxbar = impout1.xbar             # solution value"
},

{
    "location": "EAGOSolver/hp.html#",
    "page": "High Performance Builds",
    "title": "High Performance Builds",
    "category": "page",
    "text": ""
},

{
    "location": "EAGOSolver/hp.html#Configuring-EAGO-for-High-Performance-1",
    "page": "High Performance Builds",
    "title": "Configuring EAGO for High-Performance",
    "category": "section",
    "text": ""
},

{
    "location": "EAGOSolver/hp.html#Solver-Parameters-1",
    "page": "High Performance Builds",
    "title": "Solver Parameters",
    "category": "section",
    "text": "Parameter considerations for explicit optimization:Validation: The validated interval arithmetic option comes with a                 significant performance decrease but can be useful for some                 problems.\nRange Reduction: Recommend selecting an arbitrary high depth for range reduction for constrained problems.\nProbing: Recommend limiting this to the first few nodes due to the high computational cost.\nDBBT: Recommend selecting an arbitrary high depth for duality-based tightening.\nInterval Constraint Propagation: Recommended using for problems with highly nonlinear and complex constraints.Parameter considerations specifically for implicit optimization:Interval Contractor: Run roughly 5-10 interval iterations per McCormick contractor iteration. Recommend starting with 10.\nMcCormick Contractor: Limit iterations to three or fewer."
},

{
    "location": "EAGOSolver/hp.html#Lower-Bounding-Problem-1",
    "page": "High Performance Builds",
    "title": "Lower Bounding Problem",
    "category": "section",
    "text": "Currently, two (McCormick-based) lower bounding problems are available for use in both the explicit and implicit formulations: an LP solve and a solve via SNOPT. Using either problem type is recommended and the correct choice will depend on the specific problem being solve."
},

{
    "location": "EAGOSolver/hp.html#LP-Solver-Selection-1",
    "page": "High Performance Builds",
    "title": "LP Solver Selection",
    "category": "section",
    "text": "By default, EAGO uses Clp for solving linear subproblems introduced. Using a commercial linear solver is highly recommended such as Gurobi, CPLEX, or XPRESS is highly recommended. Both Gurobi and CPLEX are free for academics and installation information can be found through http://www.gurobi.com/academia/academia-center and https://www.ibm.com/developerworks/community/blogs/jfp/entry/CPLEX_Is_Free_For_Students?lang=en, respectively.  "
},

{
    "location": "EAGOSolver/hp.html#Ipopt-Build-1",
    "page": "High Performance Builds",
    "title": "Ipopt Build",
    "category": "section",
    "text": "Ipopt is the recommended solver for upper bounding problems and is supported for algorithmically optimization problems in the MathProgBase interface as well as problems defined through the JuMP interface. Ipopt\'s performance is highly dependent on the linear algebra package used (up to 30x). By default MUMPS is used. It\'s recommended that you either compile Ipopt with HSL MA57 or the Pardiso linear algebra packages with a machine specific Blas library (for Intel users the JuliaPro MKL version is recommended). For information on this, see the below links:Compiling Ipopt: https://www.coin-or.org/Ipopt/documentation/node13.html\nJulia Specifics:\nPointing Ipopt to a compiled version:\nIpopt Package Info: https://github.com/JuliaOpt/Ipopt.jl\nDiscourse discussion: https://discourse.julialang.org/t/use-ipopt-with-custom-version/9176\nIssues using Pardiso:\nUbuntu: https://github.com/JuliaOpt/Ipopt.jl/issues/106\nWindows: https://github.com/JuliaOpt/Ipopt.jl/issues/83\nHSL Website: http://www.hsl.rl.ac.uk/ipopt/\nPardiso Website: https://pardiso-project.org/"
},

{
    "location": "McCormick/overview.html#",
    "page": "Overview",
    "title": "Overview",
    "category": "page",
    "text": ""
},

{
    "location": "McCormick/overview.html#**McCormick-Operator-Capabilities**-1",
    "page": "Overview",
    "title": "McCormick Operator Capabilities",
    "category": "section",
    "text": "EAGO provides a library of McCormick relaxations in native Julia code. It supports relaxing functions using both nonsmooth McCormick relaxations (Mitsos2009), smooth McCormick relaxations (Khan2017), multi-variant McCormick relaxations (Tsoukalas2014), as well as subgradient-based interval refinement (Najman2017). For functions with arbitrarily differentiable relaxations, the differentiable can be modified by adjusting a constant value. Additionally, validated interval bounds are supported via ValidatedNumerics.jl and nonvalidated interval operators are available through use of the MCInterval{T} object included in EAGO (which is essentially just a copy of the IntervalArithmetic.jl and IntervalContractor.jl library with corrected rounding features removed)."
},

{
    "location": "McCormick/Usage.html#",
    "page": "Bounding Functions via McCormick Operators",
    "title": "Bounding Functions via McCormick Operators",
    "category": "page",
    "text": ""
},

{
    "location": "McCormick/Usage.html#**Bounding-a-function-via-smooth-McCormick-objects**-1",
    "page": "Bounding Functions via McCormick Operators",
    "title": "Bounding a function via smooth McCormick objects",
    "category": "section",
    "text": "In order to bound a function using a McCormick relaxation. You first construct structure that bounds the input variables then you construct pass these variables two a function.In the example below, convex/concave relaxations of the function f(x)=sin(2x)+exp(x)-x are calculated at x = 1 on the interval [-2,3].using EAGO\n\n# create SmoothMcCormick seed object for x = 2.0 on [1.0,3.0] for relaxing\n# a function f(x) on the interval box xIbox using mBox as a reference point\n\nf(x) = x*(x-5.0)*sin(x)\n\nx = 2.0                           # value of independent variable x\nsubx = seed_g(Float64,1,1)        # set initial subgradient of x to [1.0]\nIntv = Interval(1.0,4.0)         # define interval to relax over\n\n# create McCormick object\nSMC = SMCg{1,Interval{Float64},Float64}(x,x,subx,subx,Intv,false)\n\nfSMC = f(SMC)            # relax the function\n\ncv = fSMC.cv              # convex relaxation\ncc = fSMC.cc              # concave relaxation\ncvgrad = fSMC.cv_grad     # subgradient/gradient of convex relaxation\nccgrad = fSMC.cc_grad     # subgradient/gradient of concave relaxation\nIv = fSMC.Intv            # retrieve interval bounds of f(x) on IntvThe plotting the results we can easily generate visual the convex and concave relaxations, interval bounds, and affine bounds constructed using the subgradient at the middle of X.(Image: Figure_1)By setting the differentiability to 1, using the below command and re-plotting we arrive at the below graphset_diff_relax(1)(Image: Figure_2)This can readily be extended to multivariate functions as shown below\nset_diff_relax(0)\n\nf(x) = max(x[1],x[2])\n\nx = [2.0 1.0]                                 # values of independent variable x\nsubx = [seed_g(Float64,1,2) for i=1:2]        # set initial subgradients of x to\n                                              # [1.0, 0.0] for x[1], [0.0,1.0] for x[2]\nIntv = [Interval(-4.0,5.0),Interval(-5.0,3.0)]  # define intervals to relax over\n\n# create McCormick object\nSMC = SMCg{2,Interval{Float64},Float64}(x,x,subx,subx,Intv,false)\n\nfSMC = f(SMC)            # relax the function\n\ncv = fSMC.cv              # convex relaxation\ncc = fSMC.cc              # concave relaxation\ncvgrad = fSMC.cv_grad     # subgradient/gradient of convex relaxation\nccgrad = fSMC.cc_grad     # subgradient/gradient of concave relaxation\nIv = fSMC.Intv            # retrieve interval bounds of f(x) on Intv(Image: Figure_3)"
},

{
    "location": "McCormick/Operators.html#",
    "page": "Supported Operators",
    "title": "Supported Operators",
    "category": "page",
    "text": "For details on constructing relaxations of functions, please see the"
},

{
    "location": "McCormick/Operators.html#**Currently-supported-operators**-1",
    "page": "Supported Operators",
    "title": "Currently supported operators",
    "category": "section",
    "text": "The operators currently supported are listed below. The operators with a check box have been subject to a large degree of scrutiny and are near optimal implementations."
},

{
    "location": "McCormick/Operators.html#**Univariate-McCormick-Operators**-1",
    "page": "Supported Operators",
    "title": "Univariate McCormick Operators",
    "category": "section",
    "text": "Arbitrarily differentiable relaxations can be constructed for the following operators:[x] Inverse (inv)\n[x] Logarithms (log, log2, log10)\n[x] Exponential Functions (exp, exp2, exp10)\n[x] Square Root (sqrt)\n[x] Absolute Value (abs)Both nonsmooth and Whitney-1 (once differentiable) relaxations are supported:[x] Step Functions (step, sign)\n[x] Trignometric Functions (sin, cos, tan)\n[x] Inverse Trignometric Functions (asin, acos, atan)\n[x] Hyperbolic Functions (sinh, cosh, tanh)\n[x] Inverse Hyperbolic Functions (asinh, acosh, atanh)"
},

{
    "location": "McCormick/Operators.html#**Bivariate-Operators:-McCormick-and-McCormick**-1",
    "page": "Supported Operators",
    "title": "Bivariate Operators: McCormick & McCormick",
    "category": "section",
    "text": "The following bivariant operators are supported for two SMCg objects. Both nonsmooth and Whitney-1 (once differentiable) relaxations are supported.[x] multiplication (*)\n[x] division (/)Arbitrarily differentiable relaxations can be constructed for the following operators:[x] addition (+)\n[x] subtraction (-)\n[x] minimization (min)\n[x] maximization (max)"
},

{
    "location": "McCormick/Operators.html#**Bivariate-Operators:-McCormick-and-(Integer-or-Float)**-1",
    "page": "Supported Operators",
    "title": "Bivariate Operators: McCormick & (Integer or Float)",
    "category": "section",
    "text": "Arbitrarily differentiable relaxations can be constructed for the following operators:[x] addition (+)\n[x] subtraction (-)\n[x] multiplication (*)\n[x] division (/)\n[x] minimization (min)\n[x] maximization (max)\n[x] Exponentiation (pow, ^)"
},

{
    "location": "McCormick/Options.html#",
    "page": "Options available for McCormick Operators",
    "title": "Options available for McCormick Operators",
    "category": "page",
    "text": ""
},

{
    "location": "McCormick/Options.html#Using-fast-vs.-correctly-rounded-intervals-1",
    "page": "Options available for McCormick Operators",
    "title": "Using fast vs. correctly rounded intervals",
    "category": "section",
    "text": "Correctly rounded versus fast intervals are selected by setting the type V, in the SMCg{N,V,T} object to either:Interval{T} for correctly rounded calculations\nMCInterval{T} for fast calculations without corrected rounding"
},

{
    "location": "McCormick/Options.html#EAGO.set_diff_relax-Tuple{Integer}",
    "page": "Options available for McCormick Operators",
    "title": "EAGO.set_diff_relax",
    "category": "method",
    "text": "set_diff_relax(val::Integer)\n\nSet differentiability of relaxations used.\n\n\n\n"
},

{
    "location": "McCormick/Options.html#Setting-the-differentiability-of-the-McCormick-relaxation-1",
    "page": "Options available for McCormick Operators",
    "title": "Setting the differentiability of the McCormick relaxation",
    "category": "section",
    "text": "set_diff_relax(val::Integer)"
},

{
    "location": "McCormick/Options.html#EAGO.set_multivar_refine-Tuple{Any,Any}",
    "page": "Options available for McCormick Operators",
    "title": "EAGO.set_multivar_refine",
    "category": "method",
    "text": "set_multivar_refine(bool,tol)\n\nSet flag for using using multivariant MC relaxations and their tolerance.\n\n\n\n"
},

{
    "location": "McCormick/Options.html#Using-nonsmooth-multivariate-McCormick-relaxations-1",
    "page": "Options available for McCormick Operators",
    "title": "Using nonsmooth multivariate McCormick relaxations",
    "category": "section",
    "text": "set_multivar_refine(bool,tol)"
},

{
    "location": "McCormick/Options.html#EAGO.set_subgrad_refine-Tuple{Any}",
    "page": "Options available for McCormick Operators",
    "title": "EAGO.set_subgrad_refine",
    "category": "method",
    "text": "set_subgrad_refine(val)\n\nSet flag for using subgradient refinement of interval bounds.\n\n\n\n"
},

{
    "location": "McCormick/Options.html#Using-subgradient-refinement-1",
    "page": "Options available for McCormick Operators",
    "title": "Using subgradient refinement",
    "category": "section",
    "text": "set_subgrad_refine(val)"
},

{
    "location": "McCormick/Options.html#EAGO.set_tolerance-Tuple{Float64}",
    "page": "Options available for McCormick Operators",
    "title": "EAGO.set_tolerance",
    "category": "method",
    "text": "set_tolerance(val::Float64)\n\nSet tolerance for used in envelope calculations.\n\n\n\n"
},

{
    "location": "McCormick/Options.html#Setting-the-tolerance-used-in-envelope-calculations-1",
    "page": "Options available for McCormick Operators",
    "title": "Setting the tolerance used in envelope calculations",
    "category": "section",
    "text": "set_tolerance(val::Float64)"
},

{
    "location": "McCormick/Options.html#EAGO.set_iterations-Tuple{Integer}",
    "page": "Options available for McCormick Operators",
    "title": "EAGO.set_iterations",
    "category": "method",
    "text": "set_iterations(val::Integer)\n\nSet iterations for used in envelope calculations.\n\n\n\n"
},

{
    "location": "McCormick/Options.html#Setting-the-number-of-iterations-used-for-the-envelope-calculation-1",
    "page": "Options available for McCormick Operators",
    "title": "Setting the number of iterations used for the envelope calculation",
    "category": "section",
    "text": "set_iterations(val::Integer)"
},

{
    "location": "McCormick/back.html#",
    "page": "Back-end",
    "title": "Back-end",
    "category": "page",
    "text": ""
},

{
    "location": "McCormick/back.html#EAGO.SMCg",
    "page": "Back-end",
    "title": "EAGO.SMCg",
    "category": "type",
    "text": "SMCg{N,V,T<:AbstractFloat}\n\nSMCg is the smooth McCormick (w/ gradient) structure which is used to overload standard calculations. The fields are:\n\ncc::T: Concave relaxation\ncv::T: Convex relaxation\ncc_grad::SVector{N,T}: (Sub)gradient of concave relaxation\ncv_grad::SVector{N,T}: (Sub)gradient of convex relaxation\nIntv::V: Interval bounds\ncnst::Bool: Flag for whether the bounds are constant\n\n\n\n"
},

{
    "location": "McCormick/back.html#**The-McCormick-relaxation-type**-1",
    "page": "Back-end",
    "title": "The McCormick relaxation type",
    "category": "section",
    "text": "The McCormick relaxation library implements the smooth McCormick with imbedded (sub)gradient structure: SMCg{N,V,T}.    SMCg{N,V,T<:AbstractFloat}"
},

{
    "location": "McCormick/back.html#EAGO.HybridMC",
    "page": "Back-end",
    "title": "EAGO.HybridMC",
    "category": "type",
    "text": "HybridMC\n\nDefines the hybridMC type used for constructing nonstandard McCormick relaxations holds the SMC type.\n\n\n\n"
},

{
    "location": "McCormick/back.html#EAGO.Tighten_Subgrad",
    "page": "Back-end",
    "title": "EAGO.Tighten_Subgrad",
    "category": "function",
    "text": "Tighten_Subgrad\n\nCuts the interval bounds on the HybridMC object based on affine relaxation bounds if they are tighter.\n\n\n\n"
},

{
    "location": "McCormick/back.html#**The-Hybrid-McCormick-relaxation-type**-1",
    "page": "Back-end",
    "title": "The Hybrid McCormick relaxation type",
    "category": "section",
    "text": "The Hybrid McCormick object is used to    HybridMC    Tighten_Subgrad"
},

{
    "location": "BranchBound/Overview.html#",
    "page": "Overview",
    "title": "Overview",
    "category": "page",
    "text": "This subpart is meant to provide a flexible framework for implementing spatial branch-and-bound based optimization routines in Julia. All components of the branch-and-bound routine can be customized by the individual user: lower bounding problem, upper bounding problem. The branch and bound routine consists of a main solve algorithm that executes as depicted in the flowchart below. Routines for setting the objects to implement standard B&B routines are also provided using a set_default!() function.(Image: BnBChart1)The preprocessing routine has inputs (feas,X,UBD,k,d,opt) and outputs feas::Bool,X::Vector{Interval{Float64}}. The initial feasibility flag is feas, the bounds on the variables are X, the current upper bound is UBD, the iteration number is k, the node depth is d, and a solver option storage object is opt.\nThe lower bounding routine has inputs (X,k,d,opt,UBDg) and provides outputs (val,soln,feas,Lsto). The value of the subproblem is val, the solution of the subproblem is soln, it\'s feasibility is feas, and Lsto is a problem information storage object.\nThe upper bounding routine has inputs (X,k,d,opt,UBDg) and provides outputs (val,soln,feas,Usto). he value of the subproblem is val, the solution of the subproblem is soln, it\'s feasibility is feas, and Uto is a problem information storage object.\nThe postprocessing routine has inputs (feas,X,k,d,opt,Lsto,Usto,LBDg,UBDg) and outputs feas::Bool,X::Vector{Interval{Float64}}.\nThe repeat check has inputs (s,m,X0,X) where s::BnBSolver is a solver object, m::BnBModel is a model object, X0::Vector{Interval{Float64}} are node bounds after preprocessing, and X::Vector{Interval{Float64}} are the node bounds generated after postprocessing. Returns a boolean.\nThe bisection function has inputs (s,m,X) where s::BnBSolver is a solver object, m::BnBModel is a model object, and X::Vector{Interval{Float64}} is the box to bisect. It returns two boxes.\nThe termination check has inputs (s,m,k) where s::BnBSolver is a solver object, m::BnBModel is a model object, and k::Int64 is the iteration number. Returns a boolean.\nThe convergence check has inputs (s,UBDg,LBD) where s::BnBSolver is a solver object, UBDg is the global upper bound, and LBD is the lower bound."
},

{
    "location": "BranchBound/usage.html#",
    "page": "Usage",
    "title": "Usage",
    "category": "page",
    "text": ""
},

{
    "location": "BranchBound/usage.html#Example-1-Setup-and-Solve-a-Basic-Problem-1",
    "page": "Usage",
    "title": "Example 1 - Setup and Solve a Basic Problem",
    "category": "section",
    "text": "In the below example, we solve for minima of f(x)=x_1+x_2^2 on the domain x_1 in -11, x_2 in 29. Natural interval extensions are used to compute the upper and lower bounds. The natural interval extensions are provided by the Validated Numerics package.First, we create a BnBModel object which contains all the relevant problem info and a BnBSolver object that contains all nodes and their associated values. We specify default conditions for the Branch and Bound problem. Default conditions are a best-first search, relative width bisection, normal verbosity, a maximum of 1E6 nodes, an absolute tolerance of 1E-6, and a relative tolerance of 1E-3.\nusing EAGO\nusing ValidatedNumerics\nb = [Interval(-1,1),Interval(1,9)]\na = BnBModel(b)\nc = BnBSolver()\nEAGO.set_to_default!(c)\nc.BnB_atol = 1E-4\nNext, the lower and upper bounding problems are defined. These problems must return a tuple containing the upper/lower value, a point corresponding the upper/lower value, and the feasibility of the problem. We then set the lower/upper problem of the BnBModel object and solve the BnBModel & BnBSolver pair.\nfunction ex_LBP(X,k,pos,opt,temp)\n  ex_LBP_int = @interval X[1]+X[2]^2\n  return ex_LBP_int.lo, mid.(X), true, []\nend\nfunction ex_UBP(X,k,pos,opt,temp)\n  ex_UBP_int = @interval X[1]+X[2]^2\n  return ex_UBP_int.hi, mid.(X), true, []\nend\n\nc.Lower_Prob = ex_LBP\nc.Upper_Prob = ex_UBP\n\nouty = solveBnB!(c,a)\nThe solution is then returned in b.soln and b.UBDg is it\'s value. The corresponding output displayed to the console is given below.(Image: BnBChart2)"
},

{
    "location": "BranchBound/usage.html#Example-2-Adjust-Solver-Tolerances-1",
    "page": "Usage",
    "title": "Example 2 - Adjust Solver Tolerances",
    "category": "section",
    "text": "The absolute tolerance can be adjusted as shown below\njulia> a.BnB_tol = 1E-4\nThe relative tolerance can be changed in a similar manner\njulia> a.BnB_rtol = 1E-3\n"
},

{
    "location": "BranchBound/usage.html#Example-3-Select-Alternative-Search-Routines-1",
    "page": "Usage",
    "title": "Example 3 - Select Alternative Search Routines",
    "category": "section",
    "text": "In the above problem, the search routine could be set to a breadth-first or depth-first routine by using the set_Branch_Scheme command\njulia> EAGO.set_Branch_Scheme!(a,\"breadth\")\njulia> EAGO.set_Branch_Scheme!(a,\"depth\")\n"
},

{
    "location": "BranchBound/usage.html#Example-4-Adjust-Information-Printed-1",
    "page": "Usage",
    "title": "Example 4 - Adjust Information Printed",
    "category": "section",
    "text": "In order to print, node information in addition to iteration information the verbosity of the BnB routine can be set to full as shown below\njulia> EAGO.set_Verbosity!(a,\"Full\")\nSimilarly, if one wishes to suppress all command line outputs, the verbosity can be set to none.\njulia> EAGO.set_Verbosity!(a,\"None\")\n"
},

{
    "location": "BranchBound/API.html#",
    "page": "API",
    "title": "API",
    "category": "page",
    "text": ""
},

{
    "location": "BranchBound/API.html#EAGO.set_Branch_Scheme!-Tuple{EAGO.BnBSolver,String}",
    "page": "API",
    "title": "EAGO.set_Branch_Scheme!",
    "category": "method",
    "text": "set_Branch_Scheme!(x::BnBSolver,BM::String)\n\nSets the search scheme to \"best\", \"breadth\", or \"depth\" first schemes.\n\n\n\n"
},

{
    "location": "BranchBound/API.html#EAGO.set_Bisect_Func!-Tuple{EAGO.BnBSolver,String,Int64}",
    "page": "API",
    "title": "EAGO.set_Bisect_Func!",
    "category": "method",
    "text": "set_Bisect_Func!(x::BnBSolver,BF::String,nx::Int64)\n\nSets the bisection function to BF = \"relative midpoint\" or BF = \"absolute midpoint\" and disregards the first nx components of the interval box storage.\n\n\n\n"
},

{
    "location": "BranchBound/API.html#Using-preset-schemes-1",
    "page": "API",
    "title": "Using preset schemes",
    "category": "section",
    "text": "EAGO\'s branch and bound framework natively supports multiple common branching schemes. Namely, best-first, depth-first, and breadth-first branching schemes are all supported. The solver can be set to use any of these schemes using the set_Branch_Scheme! function below.set_Branch_Scheme!(x::BnBSolver,BM::String)Additionally, common modes of bisection are included as well. Specifically, the user can bisect the function using or absolute width bisection. For the included implicit bounding routines, the user can specify that the first nx dimensions are always ignored.set_Bisect_Func!(x::BnBSolver,BF::String,nx::Int64)"
},

{
    "location": "BranchBound/API.html#EAGO.set_Verbosity!-Tuple{EAGO.BnBSolver,String}",
    "page": "API",
    "title": "EAGO.set_Verbosity!",
    "category": "method",
    "text": "set_Verbosity!(x::BnBSolver,VB::String)\n\nSets the verbosity (console output) to either \"None\", \"Normal\", or \"Full\".\n\n\n\n"
},

{
    "location": "BranchBound/API.html#Setting-the-level-of-output-1",
    "page": "API",
    "title": "Setting the level of output",
    "category": "section",
    "text": "Currently, the branch and bound solver supports three levels of output: \"None\", \"Normal\", and \"Full\". The \"Normal\" level of output shows all iteration statistics and the final solution on termination. The \"Full\" level of output shows addition information about the \"Node\" being processed and the lower/upper bounding problems being solved.    set_Verbosity!(x::BnBSolver,VB::String)"
},

{
    "location": "BranchBound/API.html#EAGO.set_to_default!-Tuple{EAGO.BnBSolver}",
    "page": "API",
    "title": "EAGO.set_to_default!",
    "category": "method",
    "text": "set_to_default!(x::BnBSolver)\n\nReturns the B&B solver to the default settings (does not include problems or processing routines).\n\n\n\n"
},

{
    "location": "BranchBound/API.html#Returning-the-solver-to-default-settings.-1",
    "page": "API",
    "title": "Returning the solver to default settings.",
    "category": "section",
    "text": "    set_to_default!(x::BnBSolver)"
},

{
    "location": "BranchBound/API.html#EAGO.solveBnB!-Tuple{EAGO.BnBSolver,EAGO.BnBModel}",
    "page": "API",
    "title": "EAGO.solveBnB!",
    "category": "method",
    "text": "solveBnB!(x::BnBSolver,y::BnBModel)\n\nSolves the branch and bound problem with the input model and solver object.\n\n\n\n"
},

{
    "location": "BranchBound/API.html#Solving-applying-the-Branch-and-Bound-algorithm.-1",
    "page": "API",
    "title": "Solving applying the Branch and Bound algorithm.",
    "category": "section",
    "text": "    solveBnB!(x::BnBSolver,y::BnBModel)"
},

{
    "location": "BranchBound/API.html#MathProgBase.SolverInterface.getsolution-Tuple{EAGO.BnBModel}",
    "page": "API",
    "title": "MathProgBase.SolverInterface.getsolution",
    "category": "method",
    "text": "getsolution(x::BnBModel)\n\nReturns the solution stored in the BnBModel.\n\n\n\n"
},

{
    "location": "BranchBound/API.html#MathProgBase.SolverInterface.getobjval-Tuple{EAGO.BnBModel}",
    "page": "API",
    "title": "MathProgBase.SolverInterface.getobjval",
    "category": "method",
    "text": "getobjval(x::BnBModel)\n\nReturns the objective value stored in BnBModel (global upper bound).\n\n\n\n"
},

{
    "location": "BranchBound/API.html#MathProgBase.SolverInterface.getobjbound-Tuple{EAGO.BnBModel}",
    "page": "API",
    "title": "MathProgBase.SolverInterface.getobjbound",
    "category": "method",
    "text": "getobjbound(x::BnBModel)\n\nReturns the objective value stored in BnBModel (global upper bound).\n\n\n\n"
},

{
    "location": "BranchBound/API.html#EAGO.getfeasibility-Tuple{EAGO.BnBModel}",
    "page": "API",
    "title": "EAGO.getfeasibility",
    "category": "method",
    "text": "getfeasibility(x::BnBModel)\n\nReturns feasibility of problem (feasible point found?).\n\n\n\n"
},

{
    "location": "BranchBound/API.html#EAGO.LBDtime-Tuple{EAGO.BnBModel}",
    "page": "API",
    "title": "EAGO.LBDtime",
    "category": "method",
    "text": "LBDtime(x::BnBModel)\n\nReturns time spent solving lower bounding problem.\n\n\n\n"
},

{
    "location": "BranchBound/API.html#EAGO.UBDtime-Tuple{EAGO.BnBModel}",
    "page": "API",
    "title": "EAGO.UBDtime",
    "category": "method",
    "text": "UBDtime(x::BnBModel)\n\nReturns time spent solving upper bounding problem.\n\n\n\n"
},

{
    "location": "BranchBound/API.html#Accessing-info-from-solved-model-1",
    "page": "API",
    "title": "Accessing info from solved model",
    "category": "section",
    "text": "    getsolution(x::BnBModel)\n    getobjval(x::BnBModel)\n    getobjbound(x::BnBModel)\n    getfeasibility(x::BnBModel)\n    LBDtime(x::BnBModel)\n    UBDtime(x::BnBModel)"
},

{
    "location": "BranchBound/types.html#",
    "page": "Types",
    "title": "Types",
    "category": "page",
    "text": ""
},

{
    "location": "BranchBound/types.html#EAGO.BnBModel",
    "page": "Types",
    "title": "EAGO.BnBModel",
    "category": "type",
    "text": "BnBModel\n\nStores attributes of stack used to solve BnB problem. Has the following fields:\n\nInit_Box::Vector{Interval{Float64}}:        stores initial interval box used\nbox::Vector{Vector{Interval{Float64}}}      interval box storage stack\nInit_Integer::Vector{Vector{Int64}}         initial integer range\nintegers::Vector{Vector{Vector{Int64}}}     integer range storage stack\nLBD::Vector{Float64}:                       lower bounds associated with each stack item\nUBD::Vector{Float64}:                       Upper bounds associated with each stack item\nid::Vector{Int64}:                          Node ID for each stack item\npos::Vector{Int64}:                         Position in BnB Tree for each stack item\nLBDg::Float64:                              Global Lower Bound\nUBDg::Float64:                              Global Upper Bound\nLBDg_hist::Vector{Float64}:                 Value history LBD problem\nUBDg_hist::Vector{Float64}:                 Value history UBD problem\nLBDgtime::Vector{Float64}:                  Run time history LBD problem\nUBDgtime::Vector{Float64}:                  Run time history UBD problem\nPretime::Vector{Float64}:                   Run time history preprocessing\nPosttime::Vector{Float64}:                  Run time history postprocessing\nmax_id::Int64:                              Max node used\npstar::Vector{Interval{Float64}}:           IntervalBox with solution\nsoln::Vector{Float64}:                      Storage for solution\nsoln_val::Float64:                          Solution value found\nfirst_fnd::Bool:                            Has a solution been found\nfeas_fnd::Bool:                             Has a feasible point been found\nfirst_num::Int64:                           Iteration at which first solution found\nlbcnt::Int64:                               number of lower bounding problems solved\nubcnt::Int64:                               number of upper bounding problems solved\n\n\n\n"
},

{
    "location": "BranchBound/types.html#EAGO.BnBModel-Tuple{Array{IntervalArithmetic.Interval{Float64},1}}",
    "page": "Types",
    "title": "EAGO.BnBModel",
    "category": "method",
    "text": "BnBModel(X::Vector{Interval{Float64}})\n\nInitializes a BnBModel with .Init_Box = X and .box = [X].\n\n\n\n"
},

{
    "location": "BranchBound/types.html#BnBModel-1",
    "page": "Types",
    "title": "BnBModel",
    "category": "section",
    "text": "The BnBModel structure contains all information used over the course of the Branch and Bound.BnBModelAs a default, the model constructor initializes with Vector{Interval{Float64}} storage type.BnBModel(X::Vector{Interval{Float64}})"
},

{
    "location": "BranchBound/types.html#EAGO.BnBSolver",
    "page": "Types",
    "title": "EAGO.BnBSolver",
    "category": "type",
    "text": "BnBSolver\n\nStores solver specific functions used to solve BnB problem. Has the following fields:\n\nLower_Prob::Any:        Stores lower problem function (default = [])\nUpper_Prob::Any:        Stores upper problem function (default = [])\nPreprocess::Any:        Stores preprocessing function (default = [])\nPostprocess::Any:       Stores postprocessing function (default = [])\nTerm_Check::Any:        Stores termination check function (default = \'Term_Check\')\nBranch_Sto::Any:        Stores branching function (default = \'BM_depth_best!\')\nNode_Select::Any:       Stores node selection function (default = \'NS_best\')\nBisect_Func::Any:       Stores branching function (default = \'Bisect_Rel\')\nVerbosity::String:      Stores output selection (default = \"Normal\")\nmax_iter::Number:       max number of iterations (default = \"Inf\")\niter_lim::Bool:         determines if iteration limit is checked (default = false)\nmax_nodes::Int64:       max number of nodes to store in memory (default = 1E6)\nBnB_atol::Float64:      absolute tolerance for BnB (default = 1E-4)\nBnB_rtol::Float64:      relative tolerance for BnB (default = 1E-4)\nitr_intv::Int64:        number of iterations to skip between printing iteration summary (default = 20)\nhdr_intv::Int64:        number of iterations to skip between printing header (default = 1)\nconverged::Any:         convergence criterion (default = Conv_Check)\nBnB_digits::Int64:      digits displayed before decimal (default = 3)\nhist_return::Bool:      returns LBD, UBD array and time vector (default = false)\nopt::Any:               optional storage array (default = [])\nexhaust::Bool:          exhaustive search? (default = false)\ntarget_upper::Float64:  required upper bound (default = -Inf)\n\n\n\n"
},

{
    "location": "BranchBound/types.html#EAGO.BnBSolver-Tuple{}",
    "page": "Types",
    "title": "EAGO.BnBSolver",
    "category": "method",
    "text": "BnBSolver()\n\nInitializes solver with default parameters: best-first search, relative-width bisection, no iteration limit, 1E6 node limit, 1E-4 absolute and relative tolerances, no target upper bound for termination. No pre, or post processing nodes and no repetition.\n\n\n\n"
},

{
    "location": "BranchBound/types.html#BnBSolver-1",
    "page": "Types",
    "title": "BnBSolver",
    "category": "section",
    "text": "The BnBSolver options regarding how to solve the problem and routines used in it\'s solution.BnBSolverThe default initialization for BnBSolver is given below:BnBSolver()"
},

{
    "location": "BranchBound/back.html#",
    "page": "Back-end",
    "title": "Back-end",
    "category": "page",
    "text": ""
},

{
    "location": "BranchBound/back.html#Default-checks-1",
    "page": "Back-end",
    "title": "Default checks",
    "category": "section",
    "text": "Below are the default function for the Branch-and-Bound library. The EAGO solver populates these fields based on user inputs to the solver in order to deliver a valid nonconvex NLP solver.The default termination check and convergence check functions are described below:@docs\n    EAGO.Term_Check(x::BnBSolver,y::BnBModel,k_int::Int64)\n    EAGO.Conv_Check(x::BnBSolver,ubd::Float64,lbd::Float64)Currently, the default is to never repeat a node.@docs\n    EAGO.Repeat_Node_Default(x::BnBSolver,y::BnBModel{Interval{T}}, Xin::Vector{Interval{T}},Xout::Vector{Interval{T}})"
},

{
    "location": "BranchBound/back.html#Fathoming-1",
    "page": "Back-end",
    "title": "Fathoming",
    "category": "section",
    "text": "By default, nodes are fathomed on value dominance.@docs\n    EAGO.fathom!(y::BnBModel)"
},

{
    "location": "BranchBound/back.html#Pre-processing-and-post-processing-1",
    "page": "Back-end",
    "title": "Pre-processing and post-processing",
    "category": "section",
    "text": "By default, the pre-processing and post-processing functions simply return the input Interval/MCInterval type vector and the prior feasibility value."
},

{
    "location": "BranchBound/back.html#EAGO.Bisect_Abs",
    "page": "Back-end",
    "title": "EAGO.Bisect_Abs",
    "category": "function",
    "text": "EAGO.Bisect_Abs(S::BnBSolver,B::BnBModel{T},N::Vector{T})\n\nReturns two interval boxes \'X1,X2\' created by bisecting \'N\' in the highest width dimension.\n\n\n\n"
},

{
    "location": "BranchBound/back.html#EAGO.Bisect_Rel",
    "page": "Back-end",
    "title": "EAGO.Bisect_Rel",
    "category": "function",
    "text": "EAGO.Bisect_Rel(S::BnBSolver,B::BnBModel,N::Vector{Interval{Float64}})\n\nReturns two interval boxes \'X1,X2\' created by bisecting \'N\' in the highest width dimension after scaling by initial box size.\n\n\n\n"
},

{
    "location": "BranchBound/back.html#EAGO.Bisect_Abs_Imp",
    "page": "Back-end",
    "title": "EAGO.Bisect_Abs_Imp",
    "category": "function",
    "text": "EAGO.Bisect_Abs_Imp(S::BnBSolver,B::BnBModel,N::Vector{Interval{Float64}})\n\nReturns two interval boxes \'X1,X2\' created by bisecting \'N\' in the highest width dimension greater than \'nx\'.\n\n\n\n"
},

{
    "location": "BranchBound/back.html#EAGO.Bisect_Rel_Imp",
    "page": "Back-end",
    "title": "EAGO.Bisect_Rel_Imp",
    "category": "function",
    "text": "EAGO.Bisect_Rel_Imp(S::BnBSolver,B::BnBModel,N::Vector{Interval{Float64}})\n\nReturns two interval boxes \'X1,X2\' created by bisecting \'N\' in the highest width dimension greater than \'nx\' after scaling by initial box size.\n\n\n\n"
},

{
    "location": "BranchBound/back.html#Bisection-Methods-1",
    "page": "Back-end",
    "title": "Bisection Methods",
    "category": "section",
    "text": "Method for absolute width bisection on all dimension in stack:    EAGO.Bisect_AbsMethod for relative width bisection on all dimension in stack:    EAGO.Bisect_RelMethod for absolute width bisection ignore first nx dimensions:    EAGO.Bisect_Abs_ImpMethod for relative width bisection ignore first nx dimensions::    EAGO.Bisect_Rel_Imp"
},

{
    "location": "BranchBound/back.html#EAGO.BM_breadth!",
    "page": "Back-end",
    "title": "EAGO.BM_breadth!",
    "category": "function",
    "text": "EAGO.BM_breadth!\n\nTakes the following inputs: (S::BnBSolver,B::BnBModel{T},tL::Float64,tU::Float64,X1::Vector{T},X2::Vector{T},pos::Int64)\n\nStores two interval boxes X1,X2 to the bottom of the stack along with their respective lower, tL and upper bounds, tU and their position number in the BnB tree. Also, assigns node numbers.\n\n\n\n"
},

{
    "location": "BranchBound/back.html#EAGO.BM_depth_best!",
    "page": "Back-end",
    "title": "EAGO.BM_depth_best!",
    "category": "function",
    "text": "EAGO.BM_depth_best!\n\nTakes the following inputs: (S::BnBSolver,B::BnBModel,tL::Float64,tU::Float64,                              X1::Vector{Interval{Float64}},X2::Vector{Interval{Float64}},                              pos::Int64)\n\nStores two interval boxes X1,X2 to the top of the stack along with their respective lower, tL and upper bounds, tU and their position number in the BnB tree. Also, assigns node numbers.\n\n\n\n"
},

{
    "location": "BranchBound/back.html#EAGO.BM_Single!",
    "page": "Back-end",
    "title": "EAGO.BM_Single!",
    "category": "function",
    "text": "EAGO.BM_Single!\n\nTakes the following inputs: (S::BnBSolver,B::BnBModel,tL::Float64,tU::Float64, X::Vector{Interval{Float64}},pos::Int64)\n\nStores interval box X to the top of the stack along with their respective lower, tL and upper bounds, tU and their position number in the BnB tree.\n\n\n\n"
},

{
    "location": "BranchBound/back.html#Storage-Methods-for-Common-Branching-Schemes-1",
    "page": "Back-end",
    "title": "Storage Methods for Common Branching Schemes",
    "category": "section",
    "text": "Node storage method for breadth-first search:    EAGO.BM_breadth!Node storage method for depth-first or best-first search:    EAGO.BM_depth_best!Node storage method for adding a single node to the top of the stack:    EAGO.BM_Single!"
},

{
    "location": "BranchBound/back.html#EAGO.NS_best-Tuple{EAGO.BnBModel}",
    "page": "Back-end",
    "title": "EAGO.NS_best",
    "category": "method",
    "text": "EAGO.NS_best(B::BnBModel)\n\nTakes a single input B::BnBModel. Selects node with the lowest upper lower bound. Returns (IntvBox,LBD,UBD,id,pos) where Intv is the interval box, LBD is the lower bound of the node, UBD is the upper bound of the node, id is the id number of the node, and pos is the position of the node in the BnB tree.\n\n\n\n"
},

{
    "location": "BranchBound/back.html#EAGO.NS_depth_breadth-Tuple{EAGO.BnBModel}",
    "page": "Back-end",
    "title": "EAGO.NS_depth_breadth",
    "category": "method",
    "text": "EAGO.NS_depth_breadth(B::BnBModel)\n\nTakes a single input B::BnBModel. Selects node on the top of the stack. Returns (IntvBox,LBD,UBD,id,pos) where Intv is the intervalbox, LBD is the lower bound of the node, UBD is the upper bound of the node, id is the id number of the node, and pos is the position of the node in the BnB tree.\n\n\n\n"
},

{
    "location": "BranchBound/back.html#Selection-Methods-for-Common-Branching-Schemes-1",
    "page": "Back-end",
    "title": "Selection Methods for Common Branching Schemes",
    "category": "section",
    "text": "Select node for best-first search:    EAGO.NS_best(B::BnBModel)Select node for depth-first or breadth-first search:    EAGO.NS_depth_breadth(B::BnBModel)"
},

{
    "location": "BranchBound/back.html#EAGO.print_int!-Tuple{EAGO.BnBSolver,Int64,Int64,Int64,Float64,Float64,Float64,Bool,Bool}",
    "page": "Back-end",
    "title": "EAGO.print_int!",
    "category": "method",
    "text": "EAGO.print_int!(B::BnBSolver,k_int::Int64,k_nod::Int64,nid::Int64,lbdp::Float64,lbd::Float64,ubd::Float64,feasL::Bool,feasU::Bool)\n\nPrints the iteration information if the Verbosity is set to \"Normal\" or \"Full\". The header is displayed every hdr_intv, the iteration info is displayed every itr_intv\n\n\n\n"
},

{
    "location": "BranchBound/back.html#EAGO.print_results!-Tuple{EAGO.BnBSolver,Float64,Any,Bool,Bool}",
    "page": "Back-end",
    "title": "EAGO.print_results!",
    "category": "method",
    "text": "EAGO.print_results!(B::BnBSolver,sol::Float64,pnt,feas::Bool,lbd_bool::Bool)\n\nPrints the results of a single bounding problem.\n\n\n\n"
},

{
    "location": "BranchBound/back.html#EAGO.print_sol!-Tuple{EAGO.BnBSolver,EAGO.BnBModel,Int64,Int64,Float64,Float64}",
    "page": "Back-end",
    "title": "EAGO.print_sol!",
    "category": "method",
    "text": "EAGO.print_sol!(x::BnBSolver,y::BnBModel,ubdcnt::Int64,lbdcnt::Int64,ubdtime::Float64,lbdtime::Float64)\n\nPrints solution information for the B&B problem. Displays first node found, solution value, solution, and time spent solving subproblems.\n\n\n\n"
},

{
    "location": "BranchBound/back.html#Functions-for-generating-console-displayed-1",
    "page": "Back-end",
    "title": "Functions for generating console displayed",
    "category": "section",
    "text": "    EAGO.print_int!(B::BnBSolver,k_int::Int64,k_nod::Int64,nid::Int64,lbdp::Float64,lbd::Float64, ubd::Float64,feasL::Bool,feasU::Bool)    EAGO.print_results!(B::BnBSolver,sol::Float64,pnt,feas::Bool,lbd_bool::Bool)    EAGO.print_sol!(x::BnBSolver,y::BnBModel,ubdcnt::Int64,lbdcnt::Int64,ubdtime::Float64,lbdtime::Float64)"
},

{
    "location": "Dev/future.html#",
    "page": "Future Plans",
    "title": "Future Plans",
    "category": "page",
    "text": ""
},

{
    "location": "Dev/future.html#Current-Activity:-1",
    "page": "Future Plans",
    "title": "Current Activity:",
    "category": "section",
    "text": "Upgrade the solver to use MathOptInteface.\nAdd the ability to solve mixed-integer nonlinear programs via a Branch-and-Cut methodology.\nAdding support to support mixed relaxation types based on expressions and subexpressions.\nAdd support for relaxing logical-statements.\nA GUI interface for flowsheeting problems.\nSupport for propagating McCormick relaxations as expressions and taking higher order derivatives of the relaxations (port the previous EAGOSmoothMcCormick.jl package over the EAGO.jl package)"
},

{
    "location": "Dev/future.html#Other-things-on-the-wishlist-(Not-actively-being-worked-on-but-certainly-desirable):-1",
    "page": "Future Plans",
    "title": "Other things on the wishlist (Not actively being worked on but certainly desirable):",
    "category": "section",
    "text": "Implement the interval constraint propagation scheme presented in Vu 2008. Forimproved convergences.Implement dedicated single-precision linear algebra routines for intervals andMcCormick operators.Dedicated handling for linear and quadratic constraints.\nA parametric bisection routine will be updated that can divide the (X,P)space into a a series of boxes that all contain unique branches of the implicit  function p->y(p).Additional support for structured linear algebra solution routines with non-standard types.\nProvide a better interface the nonconvex semi-infinite programs solvers(JuMPeR extension?).Add additional McCormick relaxations.\nProvide means of relaxing Julia code that isn\'t flexibly typed (Cassette.jl?)"
},

{
    "location": "Dev/contributing.html#",
    "page": "How to Contribute",
    "title": "How to Contribute",
    "category": "page",
    "text": ""
},

{
    "location": "Dev/contributing.html#How-to-Contribute-1",
    "page": "How to Contribute",
    "title": "How to Contribute",
    "category": "section",
    "text": "We\'re always happy to welcome work with additional collaborators and contributors. One of the easy ways for newcomers to contribute is by adding additional relaxations objects.If you\'re interested in contributing in larger ways, please contact:Matthew Wilhelm"
},

{
    "location": "ref.html#",
    "page": "References",
    "title": "References",
    "category": "page",
    "text": ""
},

{
    "location": "ref.html#References-1",
    "page": "References",
    "title": "References",
    "category": "section",
    "text": ""
},

{
    "location": "ref.html#Branch-and-Bound-1",
    "page": "References",
    "title": "Branch and Bound",
    "category": "section",
    "text": "Floudas, Christodoulos A. Deterministic global optimization: theory, methods and applications. Vol. 37. Springer Science & Business Media, 2013.\nHorst, Reiner, and Hoang Tuy. Global optimization: Deterministic approaches. Springer Science & Business Media, 2013."
},

{
    "location": "ref.html#Parametric-Interval-Techniques-1",
    "page": "References",
    "title": "Parametric Interval Techniques",
    "category": "section",
    "text": "E. R. Hansen and G. W. Walster. Global Optimization Using Interval Analysis. Marcel Dekker, New York, second edition, 2004.\nR. Krawczyk. Newton-algorithmen zur bestimmung con nullstellen mit fehler-schranken. Computing, 4:187–201, 1969.\nR. Krawczyk. Interval iterations for including a set of solutions. Computing, 32:13–31, 1984.\nC. Miranda. Un’osservatione su un teorema di brower. Boll. Un. Mat. Ital., 3:5–7, 1940.\nA. Neumaier. Interval Methods for Systems of Equations. Cambridge University Press, Cambridge, 1990.\nR. E. Moore. A test for existence of solutions to nonlinear systems. SIAM Journal on Numerical Analysis, 14(4):611–615, 1977."
},

{
    "location": "ref.html#Domain-Reduction-1",
    "page": "References",
    "title": "Domain Reduction",
    "category": "section",
    "text": "Benhamou, F., & Older, W.J. (1997). Applying interval arithmetic to real, integer, and boolean constraints. The Journal of Logic Programming, 32, 1–24.\nCaprara, A., & Locatelli, M. (2010). Global optimization problems and domain reduction strategies. Mathematical Programming, 125, 123–137.\nGleixner, A.M., Berthold, T., Müller, B., & Weltge, S. (2016). Three enhancements for optimization-based bound tightening. ZIB Report, 15–16.\nRyoo, H.S., & Sahinidis, N.V. (1996). A branch-and-reduce approach to global optimization. Journal of Global Optimization, 8, 107–139.\nSchichl, H., & Neumaier, A. (2005). Interval analysis on directed acyclic graphs for global optimization. Journal of Global Optimization, 33, 541–562.\nTawarmalani, M., & Sahinidis, N.V. (2005). A polyhedral branch-and-cut approach to global optimization. Mathematical Programming, 103, 225–249.\nVu, X., Schichl, H., & Sam-Haroud, D. (2009). Interval propagation and search on directed acyclicgraphs for numerical constraint solving. Journal of Global Optimization, 45, 499–531."
},

{
    "location": "ref.html#Generalized-McCormick-Relaxations-1",
    "page": "References",
    "title": "Generalized McCormick Relaxations",
    "category": "section",
    "text": "Chachuat, B.: MC++: a toolkit for bounding factorable functions, v1.0. Retrieved 2 July 2014 https://projects.coin-or.org/MCpp (2014)A. Mitsos, B. Chachuat, and P. I. Barton. McCormick-based relaxations of algorithms.SIAM Journal on Optimization, 20(2):573–601, 2009.G. P. McCormick. Computability of global solutions to factorable nonconvex programs:Part I-Convex underestimating problems. Mathematical Programming, 10:147–175, 1976.G. P. McCormick. Nonlinear programming: Theory, Algorithms, and Applications. Wi-ley, New York, 1983.J. K. Scott, M. D. Stuber, and P. I. Barton. Generalized McCormick relaxations. Journalof Global Optimization, 51(4):569–606, 2011.Stuber, M.D., Scott, J.K., Barton, P.I.: Convex and concave relaxations of implicit functions. Optim.Methods Softw. 30(3), 424–460 (2015)A. Tsoukalas and A. Mitsos. Multivariate McCormick Relaxations. Journal of GlobalOptimization, 59:633–662, 2014."
},

{
    "location": "cite.html#",
    "page": "Citing EAGO",
    "title": "Citing EAGO",
    "category": "page",
    "text": ""
},

{
    "location": "cite.html#Citing-EAGO-1",
    "page": "Citing EAGO",
    "title": "Citing EAGO",
    "category": "section",
    "text": "A paper about the EAGO software package is currently under preparation. In the meantime, please feel free to cite the conference presentation below:Wilhelm, Matthew; Stuber, Matthew (October 2017) Easy Advanced Global\nOptimization (EAGO): An Open-Source Platform for Robust and Global Optimization\nin Julia. Presented at the AIChE Annual Meeting in Minneapolis, MN."
},

]}
