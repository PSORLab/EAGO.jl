# Options

## Using fast vs. correctly rounded intervals
By default, validated interval arithmetic is used via `IntervalArithmetic.jl`.

## Setting the differentiability of the McCormick relaxation
```@docs
    set_mc_differentiability!
```

## Using nonsmooth multivariate McCormick relaxations
```@docs
    set_multivar_refine!
```

## Setting the tolerance used in envelope calculations
```@docs
    set_tolerance!
```

## Setting the number of iterations used for the envelope calculation
```@docs
    set_iterations!
```
