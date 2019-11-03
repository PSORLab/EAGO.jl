# Relaxation Backend

## Quadratic Relaxations
```@docs
    EAGO.relax_convex_kernel
    EAGO.relax_nonconvex_kernel
    EAGO.relax_quadratic_gen_saf
    EAGO.relax_quadratic!
```

## Nonlinear Relaxation
```@docs
EAGO.relax_nlp!
EAGO.objective_cut_linear!
```

## Nonlinear Storage Structures
```@docs
    FunctionSetStorage
    SubexpressionSetStorage
```

## Nonlinear Evaluator
```@docs
    Evaluator
```

## Internal Functions Used by Evaluator
```@docs
    set_current_node!(x::Evaluator, n::NodeBB)
    eval_objective_lo(d::Evaluator)
    eval_constraint_cc(d::Evaluator, g::Vector{Float64}, y::Vector{Float64})
    eval_constraint_lo!(d::Evaluator, g::Vector{Float64})
    eval_constraint_hi!(d::Evaluator, g::Vector{Float64})
    eval_constraint_cc_grad(d::Evaluator, g, y)
    get_node_lower(d::FunctionSetStorage, i::Int64)
    get_node_upper(d::FunctionSetStorage, i::Int64)
    forward_reverse_pass(d::Evaluator, x::Vector{Float64})
```
