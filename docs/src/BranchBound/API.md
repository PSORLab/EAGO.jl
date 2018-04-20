## Using preset schemes
EAGO's branch and bound framework natively supports multiple common branching
schemes. Namely, best-first, depth-first, and breadth-first branching schemes
are all supported. The solver can be set to use any of these schemes using the
`set_Branch_Scheme!` function below.

```@docs
set_Branch_Scheme!(x::BnBSolver,BM::String)
```

Additionally, common modes of bisection are included as well. Specifically,
the user can bisect the function using or absolute width bisection. For the included
implicit bounding routines, the user can specify that the first nx dimensions are
always ignored.

```@docs
set_Bisect_Func!(x::BnBSolver,BF::String,nx::Int64)
```

## Setting the level of output

Currently, the branch and bound solver supports three levels of output: "None",
"Normal", and "Full". The "Normal" level of output shows all iteration statistics
and the final solution on termination. The "Full" level of output shows addition
information about the "Node" being processed and the lower/upper bounding problems
being solved.

```@docs
set_Verbosity!(x::BnBSolver,VB::String)
```

## Returning the solver to default settings.

```@docs
set_to_default!(x::BnBSolver)
```
