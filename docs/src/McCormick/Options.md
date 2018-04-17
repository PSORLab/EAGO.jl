## Using fast vs. correctly rounded intervals
Correctly rounded versus fast intervals are selected by setting the type `V`, in
the `SMCg{N,V,T}` object to either:
- `Interval{T}` for correctly rounded calculations
- `MCInterval{T}` for fast calculations without corrected rounding

## Setting the differentiability of the McCormick relaxation
```@docs
set_diff_relax(val::Integer)
```

## Using nonsmooth multivariate McCormick relaxations
```@docs
set_multivar_refine(bool,tol)
```

## Using subgradient refinement
```@docs
set_subgrad_refine(val)
```

## Setting the tolerance used in envelope calculations
```@docs
set_tolerance(val::Float64)
```

## Setting the number of iterations used for the envelope calculation
```@docs
set_iterations(val::Integer)
```
