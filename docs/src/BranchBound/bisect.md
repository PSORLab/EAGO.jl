## Bisection Methods
Method for absolute width bisection on all dimension in stack:

```@docs
Bisect_Abs(S::BnBSolver,B::BnBModel{T},N::Vector{T})
```

Method for relative width bisection on all dimension in stack:

```@docs
Bisect_Rel(S::BnBSolver,B::BnBModel,N::Vector{Interval{Float64}})
```

Method for absolute width bisection ignore first `nx` dimensions:

```@docs
Bisect_Abs_Imp(S::BnBSolver,B::BnBModel{T},N::Vector{T},nx::Q})
```

Method for relative width bisection ignore first `nx` dimensions::

```@docs
Bisect_Rel_Imp(S::BnBSolver,B::BnBModel{T},N::Vector{T},nx::Q
```
