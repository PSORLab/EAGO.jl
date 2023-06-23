# EAGO's Branch and Bound Routine

This component is meant to provide a flexible framework for implementing spatial branch-and-bound based optimization routines in Julia.
All components of the branch-and-bound routine can be customized by the individual user: lower bounding problem, upper bounding problem.

## Branch and Bound Node Storage

```@docs
    EAGO.NodeBB
```

The global optimizer structure holds all information relevant to branch-and-bound.

```@docs
    EAGO.GlobalOptimizer
```

# Customizable Subroutines

## Stack Management Subroutines

```@docs
    EAGO.branch_node!(t::ExtensionType, m::GlobalOptimizer)
    EAGO.select_branch_variable(t::ExtensionType, m::GlobalOptimizer)
    EAGO.select_branch_point(t::ExtensionType, m::GlobalOptimizer, i)
    EAGO.node_selection!(t::ExtensionType, m::GlobalOptimizer)
    EAGO.fathom!(t::ExtensionType, m::GlobalOptimizer)
    EAGO.initialize_stack!(t::ExtensionType, m::GlobalOptimizer)
    EAGO.single_storage!(t::ExtensionType, m::GlobalOptimizer)
```

## Internal Subproblem Status Codes and Subsolver Management

```@docs
    EAGO.RelaxResultStatus
    EAGO.LocalResultStatus
    EAGO.Incremental
    EAGO.SubSolvers
    EAGO.set_default_config!(t::ExtensionType, m::GlobalOptimizer)
```

## Main Subproblem and Termination Subroutines

```@docs
    EAGO.convergence_check(t::ExtensionType, m::GlobalOptimizer)
    EAGO.cut_condition(t::ExtensionType, m::GlobalOptimizer)
    EAGO.lower_problem!(t::ExtensionType, m::GlobalOptimizer{R,S,Q}) where {R,S,Q<:ExtensionType}
    EAGO.preprocess!(t::ExtensionType, m::GlobalOptimizer{R,S,Q}) where {R,S,Q<:ExtensionType}
    EAGO.postprocess!(t::ExtensionType, m::GlobalOptimizer)
    EAGO.repeat_check(t::ExtensionType, m::GlobalOptimizer)
    EAGO.termination_check(t::ExtensionType, m::GlobalOptimizer)
    EAGO.upper_problem!(t::ExtensionType, m::GlobalOptimizer)
    EAGO.parse_global!(t::ExtensionType, m::GlobalOptimizer)
    EAGO.optimize_hook!(t::ExtensionType, m::GlobalOptimizer)
```

## Internal Subroutines

```@docs
    EAGO.is_integer_subproblem(m)
    EAGO.is_integer_feasible_local(m::GlobalOptimizer, d)
│   EAGO.is_integer_feasible_relaxed(m::GlobalOptimizer)
    EAGO.interval_bound
    EAGO.lower_interval_bound
    EAGO.same_box(x::NodeBB,y::NodeBB, atol::Float64)
    EAGO.solve_local_nlp!(x::Optimizer)
    EAGO.set_dual!(x::Optimizer)
    EAGO.update_relaxed_problem_box!
    EAGO.reform_epigraph_min!(m::GlobalOptimizer)
    EAGO.label_fixed_variables!(m::GlobalOptimizer)
    EAGO.label_branch_variables!(m::GlobalOptimizer)
    EAGO.add_nonlinear!(m::GlobalOptimizer)
    EAGO.parse_classify_problem!(m::GlobalOptimizer)
    EAGO.local_problem_status!(t::MathOptInterface.TerminationStatusCode, r::MathOptInterface.ResultStatusCode)
```

## Functions for Generating Console Output

```@docs
    EAGO.print_iteration!
    EAGO.print_node!
    EAGO.print_results!
    EAGO.print_solution!
```

## Support for Log Output at Each Iteration

```@docs
    EAGO.Log
    EAGO.log_iteration!(x::GlobalOptimizer)
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
