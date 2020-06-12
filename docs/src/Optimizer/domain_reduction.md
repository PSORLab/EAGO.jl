# Domain Reduction

## Duality-Based Bound Tightening
Variable bound tightening based on the duality multipliers are supported.

```@docs
variable_dbbt!
```

## Special Forms
Bound tightening for linear forms, univariate quadratic forms, and
bivariate quadratic forms are also supported.

```@docs
EAGO.fbbt!
```

## Constraint Propagation
EAGO contains a constraint propagation architecture that supported forward and
reverse evaluation of set-valued functions on the directed acyclic graph (DAG).
The interval contractor and reverse McCormick relaxation-based contractors are
currently available.

```@docs
EAGO.set_constraint_propagation_fbbt!
```

## Optimization-Based Bound Tightening

EAGO makes use of an optimization-based bound tightening scheme using filtering
and greedy ordering as detailed in: Gleixner, A.M., Berthold, T., Müller, B.
et al. J Glob Optim (2017) 67: 731. https://doi.org/10.1007/s10898-016-0450-4.

```@docs
EAGO.obbt!
trivial_filtering!
aggressive_filtering!
bool_indx_diff
```
