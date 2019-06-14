# Basic Usage

## **Bounding a function via smooth McCormick objects**
In order to bound a function using a McCormick relaxation. You first construct
structure that bounds the input variables then you construct pass these variables
two a function.

In the example below, convex/concave relaxations of the function f(x)=sin(2x)+exp(x)-x
are calculated at x = 1 on the interval [-2,3].

```julia
using EAGO

# create MC object for x = 2.0 on [1.0,3.0] for relaxing
# a function f(x) on the interval Intv

f(x) = x*(x-5.0)*sin(x)

x = 2.0                                   # value of independent variable x
Intv = EAGO.IntervalType(1.0,4.0)         # define interval to relax over

# create McCormick object
xMC = MC{1}(x,Intv,1)

fSMC = f(xMC)            # relax the function

cv = fSMC.cv              # convex relaxation
cc = fSMC.cc              # concave relaxation
cvgrad = fSMC.cv_grad     # subgradient/gradient of convex relaxation
ccgrad = fSMC.cc_grad     # subgradient/gradient of concave relaxation
Iv = fSMC.Intv            # retrieve interval bounds of f(x) on Intv
```

The plotting the results we can easily generate visual the convex and concave
relaxations, interval bounds, and affine bounds constructed using the subgradient
at the middle of X.

![Figure_1](Figure_1.png)


By setting the differentiability to 1, using the below command and re-plotting we
arrive at the below graph
```julia
set_diff_relax!(1)
```

![Figure_2](Figure_2.png)

This can readily be extended to multivariate functions as shown below

```julia

set_mc_differentiability!(0)

f(x) = max(x[1],x[2])

x = [2.0 1.0]                                 # values of independent variable x
Intv = [EAGO.IntervalType(-4.0,5.0),EAGO.IntervalType(-5.0,3.0)]  # define intervals to relax over

# create McCormick object
xMC = [MC{2}(x[i], Intv[i], i) for i=1:2)]

fSMC = f(xSMC)            # relax the function

cv = fSMC.cv              # convex relaxation
cc = fSMC.cc              # concave relaxation
cvgrad = fSMC.cv_grad     # subgradient/gradient of convex relaxation
ccgrad = fSMC.cc_grad     # subgradient/gradient of concave relaxation
Iv = fSMC.Intv            # retrieve interval bounds of f(x) on Intv
```

![Figure_3](Figure_3.png)
