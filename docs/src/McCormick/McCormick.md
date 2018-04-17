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
Intv = Interval(-2,3)
x = EAGOSmoothMcCormick.SMC(1,1,Intv)
y = sin(2*x)+exp(x)-x
convex = y.cv # convex relaxation
concave = y.cc # concave relaxation
```
The plot obtained by computing y over all values of x is given below:

ADD PLOT HERE
