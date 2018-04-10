"""
    BnBSolver

Stores solver specific functions used to solve BnB problem. Has the following fields:
* `Lower_Prob::Any`:        Stores lower problem function (default = [])
* `Upper_Prob::Any`:        Stores upper problem function (default = [])
* `Preprocess::Any`:        Stores preprocessing function (default = [])
* `Postprocess::Any`:       Stores postprocessing function (default = [])
* `Term_Check::Any`:        Stores termination check function (default = 'Term_Check')
* `Branch_Sto::Any`:        Stores branching function (default = 'BM_depth_best!')
* `Node_Select::Any`:       Stores node selection function (default = 'NS_best')
* `Bisect_Func::Any`:       Stores branching function (default = 'Bisect_Rel')
* `Verbosity::String`:      Stores output selection (default = "Normal")
* `max_iter::Number`:       max number of iterations (default = "Inf")
* `iter_lim::Bool`:         determines if iteration limit is checked (default = false)
* `max_nodes::Int64`:       max number of nodes to store in memory (default = 1E6)
* `BnB_atol::Float64`:      absolute tolerance for BnB (default = 1E-4)
* `BnB_rtol::Float64`:      relative tolerance for BnB (default = 1E-4)
* `itr_intv::Int64`:        number of iterations to skip between printing iteration summary (default = 20)
* `hdr_intv::Int64`:        number of iterations to skip between printing header (default = 1)
* `converged::Any`:         convergence criterion (default = Conv_Check)
* `BnB_digits::Int64`:      digits displayed before decimal (default = 3)
* `hist_return::Bool`:      returns LBD, UBD array and time vector (default = false)
* `opt::Any`:               optional storage array (default = [])
* `exhaust::Bool`:          exhaustive search? (default = false)
* `target_upper::Float64`:  required upper bound (default = -Inf)
"""
mutable struct BnBSolver
  Lower_Prob::Any    # Stores lower problem function
  Upper_Prob::Any    # Stores upper problem function
  Preprocess::Any    # Preprocessing function
  Postprocess::Any   # Post processing function
  Repeat_Chk::Any
  Term_Check::Any    # Stores termination check function
  Branch_Sto::Any    # Stores branching function
  Node_Select::Any   # Stores node selection function
  Bisect_Func::Any   # Stores branching function
  Verbosity::String  # Stores output selection
  max_iter::Number   # max number of iterations
  iter_lim::Bool     # determines if iteration limit is checked
  max_nodes::Int64   # max number of nodes to store in memory
  BnB_atol::Float64   # absolute tolerance for BnB
  BnB_rtol::Float64  # relative tolerance for BnB
  itr_intv::Int64    # number of iterations to skip between printing iteration summary
  hdr_intv::Int64    # number of iterations to skip between printing header
  converged::Any     # convergence criterion
  BnB_digits::Int64  # digits displayed before decimal
  hist_return::Bool  # returns LBD, UBD array and time vector
  opt::Any           # optional storage array
  exhaust::Bool      # exhaustive search: find all solns or find first
  target_upper::Float64
end

"""
    BnBSolver()

Initializes solver with default parameters: best-first search, relative-width
bisection, no iteration limit, 1E6 node limit, 1E-4 absolute and relative
tolerances, no target upper bound for termination. No pre, or post
processing nodes and no repetition.
"""
BnBSolver() =  BnBSolver([],
                         [],
                         default_pre,
                         default_post,
                         Repeat_Node_Default,
                         Term_Check,
                         BM_depth_best!,
                         NS_best,
                         Bisect_Rel,
                         "Normal",
                         Inf,
                         false,
                         1E6,
                         1E-4,
                         1E-4,
                         1,
                         20,
                         Conv_Check,
                         3,
                         false,
                         [],
                         false,
                         -Inf)
