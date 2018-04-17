## **The McCormick relaxation structure**
The McCormick relaxation library implements the *smooth McCormick with imbedded (sub)gradient*
structure: `SMCg{N,V,T}`.

```@docs
SMCg{N,V,T<:AbstractFloat}
```

## **Bounding a function via smooth McCormick objects**
In order to bound a function using a McCormick relaxation. You first construct
structure that bounds the input variables then you construct pass these variables
two a function.

In the example below, convex/concave relaxations of the function f(x)=sin(2x)+exp(x)-x
are calculated at x = 1 on the interval [-2,3].
```julia
using EAGO
using IntervalArithmetic
using StaticArrays

a = seed_g(Float64,1,2)

# create SmoothMcCormick seed object for x1 = 2.0 on [1.0,3.0] for relaxing
# a function f(x1,x2) on the interval box xIbox using mBox as a reference point
x = 2.0
xIntv1 = Interval(1.0,3.0)
xIBox = SVector{2,Interval{Float64}}([xIntv1,xIntv1])
mBox = mid.(xIBox)
SMCgcalc = SMCg{2,Interval{Float64},Float64}(x,x,a,a,xIntv1,false,xIBox,mBox)
convex = SMCgcalc.cv # convex relaxation
concave = SMCgcalc.cc # concave relaxation
```
The plot obtained by computing y over all values of x is given below:

**ADD PLOT HERE on Github**
