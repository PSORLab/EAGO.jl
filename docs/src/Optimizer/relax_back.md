# Relaxation Backend

# Nonlinear Evaluation

EAGO makes use of a specialized tape structure for each function in order to compute valid
composite bounds and relaxations. The details of this are included below:

## Forward Operator Updates

### Main routine
```@docs
EAGO.forward_pass_kernel!
```

### Subroutines to evaluation expression
```@docs
 EAGO.forward_plus_narity!
 EAGO.forward_plus_binary!
 EAGO.forward_minus!
 EAGO.forward_multiply_binary!
 EAGO.forward_multiply_narity!
 EAGO.forward_divide!
 EAGO.forward_power!
 EAGO.forward_univariate_user!
 EAGO.forward_univariate_number!
 EAGO.forward_univariate_tiepnt_1!
 EAGO.forward_univariate_tiepnt_2!
 EAGO.forward_univariate_other!
 EAGO.forward_user_multivariate!
 EAGO.forward_get_subexpression!
```

## Reverse Operator Updates

### Main routine
```@docs
EAGO.reverse_pass!
EAGO.reverse_pass_kernel!
```

### Subroutines to evaluation expression
This list is not yet complete future documentation to
be added.
```@docs
EAGO.reverse_plus_binary!
EAGO.reverse_plus_narity!
EAGO.reverse_multiply_binary!
EAGO.reverse_multiply_narity!
```

## Other routines
```@docs
EAGO.set_value_post
EAGO.overwrite_or_intersect
```
