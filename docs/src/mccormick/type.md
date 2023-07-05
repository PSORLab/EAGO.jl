# Types

```@docs
McCormick.MC
McCormick.RelaxTag
```

## Constructors for MC

```@docs
MC{N,T}(y::Float64)
```

## Internal Utilities

```@docs
mid3
mid3v
mid_grad
dline_seg
seed_gradient
cut
secant
newton
McCormick.golden_section_it
McCormick.golden_section
```

## (Under Development) MCNoGrad

A handful of applications make use of McCormick relaxations directly without the need for subgradients. We are currently adding support for a McCormick `struct` which omits subgradient propagation in favor of return a MCNoGrad object and associated derivative information. This is currently under development and likely lacking key functionality.

```@docs
    McCormick.MCNoGrad
    MCNoGrad(y::Float64)
    MCNoGrad(y::Interval{Float64})
    MCNoGrad(cv::Float64, cc::Float64)
```
