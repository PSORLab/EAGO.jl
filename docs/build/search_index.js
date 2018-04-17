var documenterSearchIndex = {"docs": [

{
    "location": "intro.html#",
    "page": "Introduction",
    "title": "Introduction",
    "category": "page",
    "text": ""
},

{
    "location": "intro.html#EAGO-Easy-Advanced-Global-Optimization-in-Julia-1",
    "page": "Introduction",
    "title": "EAGO - Easy Advanced Global Optimization in Julia",
    "category": "section",
    "text": ""
},

{
    "location": "intro.html#Authors-1",
    "page": "Introduction",
    "title": "Authors",
    "category": "section",
    "text": "Matthew Wilhelm, Department of Chemical and Biomolecular Engineering,  University of Connecticut (UCONN)"
},

{
    "location": "intro.html#Overview-1",
    "page": "Introduction",
    "title": "Overview",
    "category": "section",
    "text": "EAGO is a global and robust optimization platform based on McCormick relaxations. It contains the first widely accessible global optimization routine based on generalized McCormick relaxations. With the exception of calls to local solvers and solve linear algebra functions all the routines are written in native Julia which allow. The solver is quite flexibly arranged so the end user can"
},

{
    "location": "intro.html#Why-McCormick-relaxations?-1",
    "page": "Introduction",
    "title": "Why McCormick relaxations?",
    "category": "section",
    "text": ""
},

{
    "location": "intro.html#Installing-EAGO-1",
    "page": "Introduction",
    "title": "Installing EAGO",
    "category": "section",
    "text": "EAGO is registered Julia package and can be installed by running:\njulia> Pkg.add(\"EAGO\")\n"
},

{
    "location": "EAGOSolver/starting.html#",
    "page": "Getting Started",
    "title": "Getting Started",
    "category": "page",
    "text": ""
},

{
    "location": "EAGOSolver/JuMP.html#",
    "page": "JuMP Interface",
    "title": "JuMP Interface",
    "category": "page",
    "text": ""
},

{
    "location": "EAGOSolver/MPB.html#",
    "page": "MathProgBase Inferace",
    "title": "MathProgBase Inferace",
    "category": "page",
    "text": ""
},

{
    "location": "EAGOSolver/SolverOpts.html#",
    "page": "Setting Solver Options",
    "title": "Setting Solver Options",
    "category": "page",
    "text": ""
},

{
    "location": "EAGOSolver/Implicit.html#",
    "page": "Implicit Function Handling",
    "title": "Implicit Function Handling",
    "category": "page",
    "text": ""
},

{
    "location": "EAGOSolver/SIP.html#",
    "page": "SIP Solver Interface",
    "title": "SIP Solver Interface",
    "category": "page",
    "text": ""
},

{
    "location": "McCormick/overview.html#",
    "page": "Overview",
    "title": "Overview",
    "category": "page",
    "text": ""
},

{
    "location": "McCormick/McCormick.html#",
    "page": "Bounding Functions via McCormick Operators",
    "title": "Bounding Functions via McCormick Operators",
    "category": "page",
    "text": ""
},

{
    "location": "McCormick/Options.html#",
    "page": "Options available for McCormick Operators",
    "title": "Options available for McCormick Operators",
    "category": "page",
    "text": ""
},

{
    "location": "ParamInterval/contractor.html#",
    "page": "Parametric Contractor",
    "title": "Parametric Contractor",
    "category": "page",
    "text": ""
},

{
    "location": "ParamInterval/contractor.html#Parametric-Interval-Contractor-1",
    "page": "Parametric Contractor",
    "title": "Parametric Interval Contractor",
    "category": "section",
    "text": "Provides methods for performing parametric interval calculations such as (Parametric Interval Newton/Krawczyk) as well as a series of tests to verify the (non)existence of unique enclosed functions. "
},

{
    "location": "ParamInterval/test.html#",
    "page": "Parametric Tests",
    "title": "Parametric Tests",
    "category": "page",
    "text": ""
},

{
    "location": "DomainReduction/overview.html#",
    "page": "Overview",
    "title": "Overview",
    "category": "page",
    "text": ""
},

{
    "location": "DomainReduction/RR.html#",
    "page": "Range Reduction",
    "title": "Range Reduction",
    "category": "page",
    "text": ""
},

{
    "location": "DomainReduction/DBBT.html#",
    "page": "Duality Based Bound Tightening",
    "title": "Duality Based Bound Tightening",
    "category": "page",
    "text": ""
},

{
    "location": "DomainReduction/probing.html#",
    "page": "Probing",
    "title": "Probing",
    "category": "page",
    "text": ""
},

{
    "location": "BranchBound/overview.html#",
    "page": "Overview",
    "title": "Overview",
    "category": "page",
    "text": "This package is meant to provide a flexible framework for implementing branch-and-bound based optimization routines in Julia. All components of the branch-and-bound routine can be customized by the individual user: lower bounding problem, upper bounding problem. The branch and bound routine consists of a main solve algorithm that executes as depicted in the flowchart below. Routines for setting the objects to implement standard B&B routines are also provided using a set_default!() function.ADD CHART HERE ON GITHUB!!!!!!!The preprocessing routine has inputs (feas,X,UBD,k,d,opt) and outputs feas::Bool,X::Vector{Interval{Float64}}. The initial feasibility flag is feas, the bounds on the variables are X, the current upper bound is UBD, the iteration number is k, the node depth is d, and a solver option storage object is opt.\nThe lower bounding routine has inputs (X,k,d,opt,UBDg) and provides outputs (val,soln,feas,Lsto). The value of the subproblem is val, the solution of the subproblem is soln, it\'s feasibility is feas, and Lsto is a problem information storage object.\nThe upper bounding routine has inputs (X,k,d,opt,UBDg) and provides outputs (val,soln,feas,Usto). he value of the subproblem is val, the solution of the subproblem is soln, it\'s feasibility is feas, and Uto is a problem information storage object.\nThe postprocessing routine has inputs (feas,X,k,d,opt,Lsto,Usto,LBDg,UBDg) and outputs feas::Bool,X::Vector{Interval{Float64}}.\nThe repeat check has inputs (s,m,X0,X) where s::BnBSolver is a solver object, m::BnBModel is a model object, X0::Vector{Interval{Float64}} are node bounds after preprocessing, and X::Vector{Interval{Float64}} are the node bounds generated after postprocessing. Returns a boolean.\nThe bisection function has inputs (s,m,X) where s::BnBSolver is a solver object, m::BnBModel is a model object, and X::Vector{Interval{Float64}} is the box to bisect. It returns two boxes.\nThe termination check has inputs (s,m,k) where s::BnBSolver is a solver object, m::BnBModel is a model object, and k::Int64 is the iteration number. Returns a boolean.\nThe convergence check has inputs (s,UBDg,LBD) where s::BnBSolver is a solver object, UBDg is the global upper bound, and LBD is the lower bound."
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
    "text": "In the below example, we solve for minima of f(x)=x<sub>1</sub>+x<sub>2</sub><sup>2</sup> on the domain [-1,1] by [2,9]. Natural interval extensions are used to compute the upper and lower bounds. The natural interval extensions are provided by the Validated Numerics package.First, we create a BnBModel object which contains all the relevant problem info and a BnBSolver object that contains all nodes and their associated values. We specify default conditions for the Branch and Bound problem. Default conditions are a best-first search, relative width bisection, normal verbosity, a maximum of 1E6 nodes, an absolute tolerance of 1E-6, and a relative tolerance of 1E-3.\nusing EAGO\nusing ValidatedNumerics\nb = [Interval(-1,1),Interval(1,9)]\na = BnBModel(b)\nc = BnBSolver()\nEAGO.set_to_default!(c)\nc.BnB_atol = 1E-4\nNext, the lower and upper bounding problems are defined. These problems must return a tuple containing the upper/lower value, a point corresponding the upper/lower value, and the feasibility of the problem. We then set the lower/upper problem of the BnBModel object and solve the BnBModel & BnBSolver pair.\nfunction ex_LBP(X,k,pos,opt,temp)\n  ex_LBP_int = @interval X[1]+X[2]^2\n  return ex_LBP_int.lo, mid.(X), true, []\nend\nfunction ex_UBP(X,k,pos,opt,temp)\n  ex_UBP_int = @interval X[1]+X[2]^2\n  return ex_UBP_int.hi, mid.(X), true, []\nend\n\nc.Lower_Prob = ex_LBP\nc.Upper_Prob = ex_UBP\n\nouty = solveBnB!(c,a)\nThe solution is then returned in b.soln and b.UBDg is it\'s value. The corresponding output displayed to the console is given below.ADD CHART HERE ON GITHUB!!!!!!!"
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
    "location": "BranchBound/types.html#EAGO.BnBModel",
    "page": "Types",
    "title": "EAGO.BnBModel",
    "category": "type",
    "text": "BnBModel\n\nStores attributes of stack used to solve BnB problem. Has the following fields:\n\nInit_Box::Vector{Interval{Float64}}:        stores initial interval box used\nbox::Vector{Vector{Interval{Float64}}}      interval box storage stack\nInit_Integer::Vector{Vector{Int64}}         initial integer range\nintegers::Vector{Vector{Vector{Int64}}}     integer range storage stack\nLBD::Vector{Float64}:                       lower bounds associated with each stack item\nUBD::Vector{Float64}:                       Upper bounds associated with each stack item\nid::Vector{Int64}:                          Node ID for each stack item\npos::Vector{Int64}:                         Position in BnB Tree for each stack item\nLBDg::Float64:                              Global Lower Bound\nUBDg::Float64:                              Global Upper Bound\nLBDg_hist::Vector{Float64}:                 Value history LBD problem\nUBDg_hist::Vector{Float64}:                 Value history UBD problem\nLBDgtime::Vector{Float64}:                  Run time history LBD problem\nUBDgtime::Vector{Float64}:                  Run time history UBD problem\nPretime::Vector{Float64}:                   Run time history preprocessing\nPosttime::Vector{Float64}:                  Run time history postprocessing\nmax_id::Int64:                              Max node used\npstar::Vector{Interval{Float64}}:           IntervalBox with solution\nsoln::Vector{Float64}:                      Storage for solution\nsoln_val::Float64:                          Solution value found\nfirst_fnd::Bool:                            Has a solution been found\nfeas_fnd::Bool:                             Has a feasible point been found\nfirst_num::Int64:                           Iteration at which first solution found\nlbcnt::Int64:                               number of lower bounding problems solved\nubcnt::Int64:                               number of upper bounding problems solved\n\n\n\n"
},

{
    "location": "BranchBound/types.html#EAGO.BnBSolver",
    "page": "Types",
    "title": "EAGO.BnBSolver",
    "category": "type",
    "text": "BnBSolver\n\nStores solver specific functions used to solve BnB problem. Has the following fields:\n\nLower_Prob::Any:        Stores lower problem function (default = [])\nUpper_Prob::Any:        Stores upper problem function (default = [])\nPreprocess::Any:        Stores preprocessing function (default = [])\nPostprocess::Any:       Stores postprocessing function (default = [])\nTerm_Check::Any:        Stores termination check function (default = \'Term_Check\')\nBranch_Sto::Any:        Stores branching function (default = \'BM_depth_best!\')\nNode_Select::Any:       Stores node selection function (default = \'NS_best\')\nBisect_Func::Any:       Stores branching function (default = \'Bisect_Rel\')\nVerbosity::String:      Stores output selection (default = \"Normal\")\nmax_iter::Number:       max number of iterations (default = \"Inf\")\niter_lim::Bool:         determines if iteration limit is checked (default = false)\nmax_nodes::Int64:       max number of nodes to store in memory (default = 1E6)\nBnB_atol::Float64:      absolute tolerance for BnB (default = 1E-4)\nBnB_rtol::Float64:      relative tolerance for BnB (default = 1E-4)\nitr_intv::Int64:        number of iterations to skip between printing iteration summary (default = 20)\nhdr_intv::Int64:        number of iterations to skip between printing header (default = 1)\nconverged::Any:         convergence criterion (default = Conv_Check)\nBnB_digits::Int64:      digits displayed before decimal (default = 3)\nhist_return::Bool:      returns LBD, UBD array and time vector (default = false)\nopt::Any:               optional storage array (default = [])\nexhaust::Bool:          exhaustive search? (default = false)\ntarget_upper::Float64:  required upper bound (default = -Inf)\n\n\n\n"
},

{
    "location": "BranchBound/types.html#",
    "page": "Types",
    "title": "Types",
    "category": "page",
    "text": "The BnBModel structure contains all information used over the course of the Branch and Bound.BnBModelAs a default, the model constructor initializes with Vector{Interval{Float64}} storage type.The BnBSolver options regarding how to solve the problem and routines used in it\'s solution.BnBSolver"
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
    "text": "Currently, the branch and bound solver supports three levels of output: \"None\", \"Normal\", and \"Full\". The \"Normal\" level of output shows all iteration statistics and the final solution on termination. The \"Full\" level of output shows addition information about the \"Node\" being processed and the lower/upper bounding problems being solved.set_Verbosity!(x::BnBSolver,VB::String)"
},

{
    "location": "BranchBound/API.html#EAGO.set_to_default!-Tuple{EAGO.BnBSolver}",
    "page": "API",
    "title": "EAGO.set_to_default!",
    "category": "method",
    "text": "set_to_default!(x::BnBSolver)\n\nReturns the B&B solver to the default settings (does not include problems or processing routines.\n\n\n\n"
},

{
    "location": "BranchBound/API.html#Returning-the-solver-to-default-settings.-1",
    "page": "API",
    "title": "Returning the solver to default settings.",
    "category": "section",
    "text": "set_to_default!(x::BnBSolver)"
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
    "text": "A paper about the EAGO software package is currently under preparation. In the meantime, please feel free to cite the conference presentation below:\nWilhelm, Matthew; Stuber, Matthew (October 2017) Easy Advanced Global\nOptimization (EAGO): An Open-Source Platform for Robust and Global Optimization\nin Julia. Presented at the AIChE Annual Meeting in Minneapolis, MN.\n"
},

]}
