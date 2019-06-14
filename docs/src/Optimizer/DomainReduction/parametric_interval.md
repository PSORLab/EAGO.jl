# Parametric Interval methods

For those of you unfamiliar with basic interval arithmetic, please check out the
[ValidatedNumerics.jl](https://github.com/JuliaIntervals/ValidatedNumerics.jl) package.

Parametric interval contractor methods are quite similar to their standard interval
analogs. The key exception is that interval contractor methods, e.g. Newton/Krawczyk,
bound roots of a function, the parametric methods bound functions parametrically defined by these roots.

To elaborate further, the interval contractor method provides bounds, X, on $x \in X$ s.t. $f(x) = 0$.
The parameter interval contractor methods takes a function $y = x(p)$ defined implicitly by $h(y,p) = 0$ on (X,P)
and then subsequent computes a smaller box containing a solution branch of $h(y,p) = 0$.

## Parametric Interval Contractor

EAGO includes a series of interval parametric interval contractor routines along with.
Provides methods for performing parametric interval calculations such as (Parametric Interval Newton/Krawczyk) as well as a series of tests to verify the (non)existence of unique enclosed functions.

```@docs
param_intv_contractor
```

The options used for interval contractor routines are specified by an option structure.

```@docs
parametric_interval_params
```

### Usage

An example of using the interval contractors is given below:

```julia
opt1 = parametric_interval_params(:Dense,:Krawczyk,1E-30,1E-6,2,0,100)
opt2 = parametric_interval_params(:Dense,:Newton,1E-30,1E-6,2,0,100)


P1 = [Interval(5.0,7.0)Interval(5.0,7.0)]     # define control variable bounds
Z1 = [Interval(-1.5, 0.0),Interval(0.0, 0.5)]  # define state variable bounds
Z1c = copy(Z1)

# Defines the system to contract/test
h1(z,p) = [z[1]^2+z[2]^2+p[1]*z[1]+4.0;
           z[1]+p[2]*z[2]]

# Define the jacobian
hj1(z,p) = [(z[1]+p[1]) (z[2]);
              one(p[1])  p[2]]

# Set initial flags to indicate results are unknown
Eflag = false
Iflag = false
eDflag = false

# Run the contractors
output1 = param_intv_contractor(h1,hj1,Z1,P1,Eflag,Iflag,eDflag,opt1)
output2 = param_intv_contractor(h1,hj1,Z1c,P1,Eflag,Iflag,eDflag,opt2)
```
