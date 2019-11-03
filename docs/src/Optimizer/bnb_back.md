# EAGO's Branch and Bound Routine

This component is meant to provide a flexible framework for implementing spatial branch-and-bound based optimization routines in Julia.
All components of the branch-and-bound routine can be customized by the individual user: lower bounding problem, upper bounding problem.

## Branch and Bound Node Storage
```@docs
    EAGO.NodeBB
```

## Customizable subroutines

```@docs
    EAGO.add_cut!(t::ExtensionType, x::Optimizer)
    EAGO.branch_node!(t::ExtensionType, x::Optimizer)
    EAGO.convergence_check(t::ExtensionType, x::Optimizer)
    EAGO.cut_condition(t::ExtensionType, x::Optimizer)
    EAGO.fathom!(t::ExtensionType, d::Optimizer)
    EAGO.lower_problem!(t::ExtensionType, x::Optimizer)
    EAGO.relax_objective!(t::ExtensionType, x::Optimizer, x0::Vector{Float64})
    EAGO.relax_problem!(t::ExtensionType, x::Optimizer, v::Vector{Float64}, q::Int64)
    EAGO.node_selection!(t::ExtensionType, x::Optimizer)
    EAGO.optimize_hook!(t::ExtensionType, x::Optimizer)
    EAGO.postprocess!(t::ExtensionType, x::Optimizer)
    EAGO.preprocess!(t::ExtensionType, x::Optimizer)
    EAGO.repeat_check(t::ExtensionType, x::Optimizer)
    EAGO.single_storage!(t::ExtensionType, x::Optimizer)
    EAGO.termination_check(t::ExtensionType, x::Optimizer)
    EAGO.upper_problem!(t::ExtensionType, x::Optimizer)
```

## Internal Subroutines
```@docs
    EAGO.cut_update(x::Optimizer)
    EAGO.default_nlp_heurestic(x::Optimizer, y::NodeBB)
    EAGO.interval_bound
    EAGO.interval_lower_bound!(x::Optimizer, y::NodeBB)
    EAGO.is_globally_optimal(t::MOI.TerminationStatusCode, r::MOI.ResultStatusCode)
    EAGO.is_feasible_solution(t::MOI.TerminationStatusCode, r::MOI.ResultStatusCode)
    EAGO.log_iteration!(x::Optimizer)
    EAGO.same_box(x::NodeBB,y::NodeBB, atol::Float64)
    EAGO.solve_local_nlp!(x::Optimizer)
    EAGO.set_dual!(x::Optimizer)
    EAGO.update_relaxed_problem_box!(x::Optimizer, y::NodeBB)
```

## Functions for generating console output
```@docs
    EAGO.print_iteration!
    EAGO.print_node!
    EAGO.print_results!
    EAGO.print_solution!
```
