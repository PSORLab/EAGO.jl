"""
    BnBModel
Stores attributes of stack used to solve BnB problem. Has the following fields:
* `Init_Box::Vector{Interval{Float64}}`:        stores initial interval box used
* `box::Vector{Vector{Interval{Float64}}}`      interval box storage stack
* `Init_Integer::Vector{Vector{Int64}}`         initial integer range
* `integers::Vector{Vector{Vector{Int64}}}`     integer range storage stack
* `LBD::Vector{Float64}`:                       lower bounds associated with each stack item
* `UBD::Vector{Float64}`:                       Upper bounds associated with each stack item
* `id::Vector{Int64}`:                          Node ID for each stack item
* `pos::Vector{Int64}`:                         Position in BnB Tree for each stack item
* `LBDg::Float64`:                              Global Lower Bound
* `UBDg::Float64`:                              Global Upper Bound
* `LBDg_hist::Vector{Float64}`:                 Value history LBD problem
* `UBDg_hist::Vector{Float64}`:                 Value history UBD problem
* `LBDgtime::Vector{Float64}`:                  Run time history LBD problem
* `UBDgtime::Vector{Float64}`:                  Run time history UBD problem
* `Pretime::Vector{Float64}`:                   Run time history preprocessing
* `Posttime::Vector{Float64}`:                  Run time history postprocessing
* `max_id::Int64`:                              Max node used
* `pstar::Vector{Interval{Float64}}`:           IntervalBox with solution
* `soln::Vector{Float64}`:                      Storage for solution
* `soln_val::Float64`:                          Solution value found
* `first_fnd::Bool`:                            Has a solution been found
* `feas_fnd::Bool`:                             Has a feasible point been found
* `first_num::Int64`:                           Iteration at which first solution found
* `lbcnt::Int64`:                               number of lower bounding problems solved
* `ubcnt::Int64`:                               number of upper bounding problems solved
"""
mutable struct BnBModel{V}
  Init_Box::Vector{V}
  box::Vector{Vector{V}} # interval box stack
  Init_Integer::Vector{Vector{Int64}}
  integers::Vector{Vector{Vector{Int64}}}
  LBD::Vector{Float64} # lower bounds associated with each stack
  UBD::Vector{Float64} # Upper bounds associated with each stack
  id::Vector{Int64}  # Node ID
  pos::Vector{Int64} # Position in BnB Tree
  LBDg::Float64 # Global Lower Bound
  UBDg::Float64 # Global Upper Bound
  LBDg_hist::Vector{Float64} # Value history LBD problem
  UBDg_hist::Vector{Float64} # Value history UBD problem
  LBDgtime::Vector{Float64} # Run time history LBD problem
  UBDgtime::Vector{Float64} # Run time history UBD problem
  Pretime::Vector{Float64} # Run time history LBD problem
  Posttime::Vector{Float64} # Run time history UBD problem
  max_id::Int64             # Max node used
  pstar::Vector{V}
  soln::Vector{Float64}
  soln_val::Float64
  first_fnd::Bool
  feas_fnd::Bool
  first_num::Int64
  lbcnt::Int64
  ubcnt::Int64
  lastgap::Float64
end

"""
    BnBModel(X::Vector{Interval{Float64}})

Initializes a `BnBModel` with `.Init_Box` = `X` and `.box` = `[X]`.
"""
BnBModel(X::Vector{V}) where {V} = BnBModel{V}(deepcopy(X),
                                            [deepcopy(X)],
                                            [[Int64(1)]],
                                            [[[Int64(1)]]],
                                            [-Inf],
                                            [Inf],
                                            [Int64(1)],
                                            [Int64(1)],
                                            -Inf,
                                            Inf,
                                            [-Inf],
                                            [Inf],
                                            [0.0],
                                            [0.0],
                                            [0.0],
                                            [0,0],
                                            Int64(1),
                                            deepcopy(X),
                                            [0.0],
                                            Inf,
                                            false,
                                            false,
                                            -Int64(1),
                                            Int64(0),
                                            Int64(0),
                                            Inf)
#BnBModel() = BnBModel{MCInterval{Float64}}([MCInterval(0.0,1.0)])
