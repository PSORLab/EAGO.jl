# Basic Usage

## Bounding a Function via McCormick Operators

In order to bound a function using a McCormick relaxation, you first construct a
McCormick object (`x::MC`) that bounds the input variables, and then you pass these
variables to the desired function.

In the example below, convex/concave relaxations of the function `f(x) = x(x-5)sin(x)`
are calculated at `x = 2` on the interval `[1,4]`.

```julia
using McCormick

# create MC object for x = 2.0 on [1.0,4.0] for relaxing
# a function f(x) on the interval Intv

f(x) = x*(x-5.0)*sin(x)

x = 2.0                          # value of independent variable x
Intv = Interval(1.0,4.0)         # define interval to relax over
                                 # Note that McCormick.jl reexports IntervalArithmetic.jl
                                 # and StaticArrays. So no using statement for these is
                                 # necessary.
# create McCormick object
xMC = MC{1,NS}(x,Intv,1)

fMC = f(xMC)             # relax the function

cv = fMC.cv              # convex relaxation
cc = fMC.cc              # concave relaxation
cvgrad = fMC.cv_grad     # subgradient/gradient of convex relaxation
ccgrad = fMC.cc_grad     # subgradient/gradient of concave relaxation
Iv = fMC.Intv            # retrieve interval bounds of f(x) on Intv
```

By plotting the results we can easily visualize the convex and concave
relaxations, interval bounds, and affine bounds constructed using the subgradient
at the middle of X.

![Figure_1](Figure_1.png)

If we instead use the constructor `xMC = MC{1,Diff}(x,Intv,1)` in the above code and re-plot, 
we arrive at the following graph. Note that these relaxations are differentiable, but not as
tight as the nonsmooth relaxations.

![Figure_2](Figure_2.png)

This can readily be extended to multivariate functions, for example, `f(x,y) = (4 - 2.1x^2 + (x^4)/6)x^2 + xy + (-4 + 4y^2)y^2`:

```julia
using McCormick

# initialize function
f(x,y) = (4.0 - 2.1*x^2 + (x^4)/6.0)*x^2 + x*y + (-4.0 + 4.0*y^2)*y^2

# intervals for independent variables
n = 30
X = Interval{Float64}(-2,0)
Y = Interval{Float64}(-0.5,0.5)
xrange = range(X.lo,stop=X.hi,length=n)
yrange = range(Y.lo,stop=Y.hi,length=n)

# differentiable McCormick relaxation
for (i,x) in enumerate(xrange)
    for (j,y) in enumerate(yrange)
        z = f(x,y)                  # calculate function values
        xMC = MC{1,Diff}(x,X,1)     # differentiable relaxation for x
        yMC = MC{1,Diff}(y,Y,2)     # differentiable relaxation for y
        fMC = f(xMC,yMC)            # relax the function
        cv = fMC.cv                 # convex relaxation
        cc = fMC.cc                 # concave relaxation
    end
end
```

![Figure_4](Figure_4.png)
