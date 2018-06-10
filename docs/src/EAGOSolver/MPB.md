## Using the MathProgBase Interface

The MathProgBase interface is recommended for optimizing nonsmooth objects and
problems that can't be easily represented in JuMP. A EAGO model object should
be first created from as solver

```julia
f(x) = (x[1]-5)^2 + (x[2]-3)^2
g(x) = [x[1] - x[2]]
m1 = MathProgBase.NonlinearModel(EAGO_NLPSolver())
MathProgBase.loadproblem!(m1, 2, 1, [0.0, 0.0], [10.0, 10.0],
            [0.0], [1.0], :Min, f, g)
MathProgBase.optimize!(m1)
```

## EAGO Model Object

```@docs
EAGO_NLP_Model
```

```@docs
EAGO_Inner_NLP
```

## Setup API for MathProgBase Interface


```@docs
    MathProgBase.optimize!(s::EAGO_NLP_Model)
```

## Access and Manipulation functions for MathProgBase Interface

Standard MathProgBase access functions (e.g. MathProgBase.getsolution) are extended by EAGO.
