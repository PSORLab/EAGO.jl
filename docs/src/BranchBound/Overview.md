
This subpart is meant to provide a flexible framework for implementing spatial branch-and-bound based optimization routines in Julia.
All components of the branch-and-bound routine can be customized by the individual user: lower bounding problem, upper bounding problem.
The branch and bound routine consists of a main solve algorithm that executes as depicted in the flowchart below.
Routines for setting the objects to implement standard B&B routines are also provided using a `set_default!()` function.

![BnBChart1](BnBChart1.jpg)

- The preprocessing routine has inputs `(feas,X,UBD,k,d,opt)` and outputs `feas::Bool,X::Vector{Interval{Float64}}`. The initial feasibility flag is `feas`, the bounds on the variables are `X`, the current upper bound is `UBD`, the iteration number is `k`, the node depth is `d`, and a solver option storage object is `opt`.
- The lower bounding routine has inputs `(X,k,d,opt,UBDg)` and provides outputs `(val,soln,feas,Lsto)`. The value of the subproblem is `val`, the solution of the subproblem is `soln`, it's feasibility is `feas`, and `Lsto` is a problem information storage object.
- The upper bounding routine has inputs `(X,k,d,opt,UBDg)` and provides outputs `(val,soln,feas,Usto)`. he value of the subproblem is `val`, the solution of the subproblem is `soln`, it's feasibility is `feas`, and `Uto` is a problem information storage object.
- The postprocessing routine has inputs `(feas,X,k,d,opt,Lsto,Usto,LBDg,UBDg)` and outputs `feas::Bool,X::Vector{Interval{Float64}}`.
- The repeat check has inputs `(s,m,X0,X)` where `s::BnBSolver` is a solver object, `m::BnBModel` is a model object, `X0::Vector{Interval{Float64}}` are node bounds after preprocessing, and `X::Vector{Interval{Float64}}` are the node bounds generated after postprocessing. Returns a boolean.
- The bisection function has inputs `(s,m,X)` where `s::BnBSolver` is a solver object, `m::BnBModel` is a model object, and `X::Vector{Interval{Float64}}` is the box to bisect. It returns two boxes.
- The termination check has inputs `(s,m,k)` where `s::BnBSolver` is a solver object, `m::BnBModel` is a model object, and `k::Int64` is the iteration number. Returns a boolean.
- The convergence check has inputs `(s,UBDg,LBD)` where `s::BnBSolver` is a solver object, `UBDg` is the global upper bound, and `LBD` is the lower bound.
