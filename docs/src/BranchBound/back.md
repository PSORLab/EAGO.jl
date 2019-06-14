# Backend


## Branch and Bound Node Storage
```@docs
    EAGO.NodeBB
```

## Default checks

Below are the default function for the Branch-and-Bound library. The EAGO solver
populates these fields based on user inputs to the solver in order to deliver a
valid nonconvex NLP solver.

The default repetition, termination, and convergence check functions are described below:
```@docs
    EAGO.default_repeat_check
    EAGO.default_termination_check
    EAGO.default_convergence_check
```

## Fathoming

By default, nodes are fathomed on value dominance.
```@docs
    EAGO.fathom!
```

## Bisection Methods
Method for relative width bisection on all dimension in stack:

```@docs
    EAGO.continuous_relative_bisect
```

Method for relative width bisection ignoring the last `nx` dimensions:

```@docs
    EAGO.implicit_bisection
```

## Selection Methods for Common Branching Schemes

Select node for best-first search:

```@docs
    EAGO.node_select_best!
```

## Functions for generating console displayed
```@docs
    EAGO.print_node!
    EAGO.print_iteration!
    EAGO.print_solution!
    EAGO.print_results!
```
