# EAGO's Branch and Bound Routine

This component is meant to provide a flexible framework for implementing spatial branch-and-bound based optimization routines in Julia.
All components of the branch-and-bound routine can be customized by the individual user: lower bounding problem, upper bounding problem.

## Branch and Bound Node Storage
```@docs
    EAGO.NodeBB
```

## Customizable subroutines

```@docs
    EAGO.branch_node!(t::ExtensionType, x::GlobalOptimizer{R,S,Q}) where {R,S,Q<:ExtensionType}
    EAGO.convergence_check(t::ExtensionType, x::GlobalOptimizer{R,S,Q}) where {R,S,Q<:ExtensionType}
    EAGO.fathom!(t::ExtensionType, d::GlobalOptimizer{R,S,Q}) where {R,S,Q<:ExtensionType}
    EAGO.lower_problem!(t::ExtensionType, x::GlobalOptimizer{R,S,Q}) where {R,S,Q<:ExtensionType}
    EAGO.node_selection!(t::ExtensionType, x::GlobalOptimizer{R,S,Q}) where {R,S,Q<:ExtensionType}
    EAGO.postprocess!(t::ExtensionType, x::GlobalOptimizer{R,S,Q}) where {R,S,Q<:ExtensionType}
    EAGO.preprocess!(t::ExtensionType, x::GlobalOptimizer{R,S,Q}) where {R,S,Q<:ExtensionType}
    EAGO.repeat_check(t::ExtensionType, x::GlobalOptimizer{R,S,Q}) where {R,S,Q<:ExtensionType}
    EAGO.single_storage!(t::ExtensionType, x::GlobalOptimizer{R,S,Q}) where {R,S,Q<:ExtensionType}
    EAGO.termination_check(t::ExtensionType, x::GlobalOptimizer{R,S,Q}) where {R,S,Q<:ExtensionType}
    EAGO.upper_problem!(t::ExtensionType, x::GlobalOptimizer{R,S,Q}) where {R,S,Q<:ExtensionType}
```

## Internal Subroutines
```@docs
    EAGO.interval_bound
    EAGO.lower_interval_bound
    EAGO.log_iteration!(x::Optimizer)
    EAGO.same_box(x::NodeBB,y::NodeBB, atol::Float64)
    EAGO.solve_local_nlp!(x::Optimizer)
    EAGO.set_dual!(x::Optimizer)
    EAGO.update_relaxed_problem_box!
```

## Functions for generating console output
```@docs
    EAGO.print_iteration!
    EAGO.print_node!
    EAGO.print_results!
    EAGO.print_solution!
```

## Interval Representations of Expressions
```@docs
    EAGO.AbstractEAGOConstraint
    EAGO.AffineFunctionEq
    EAGO.AffineFunctionIneq
    EAGO.BufferedQuadraticIneq
    EAGO.BufferedQuadraticEq
    EAGO.NonlinearExpression
    EAGO.BufferedNonlinearFunction
```
