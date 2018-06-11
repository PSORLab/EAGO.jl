Below are the default function for the Branch-and-Bound library. The EAGO solver
populates these fields based on user inputs to the solver in order to deliver a
valid nonconvex NLP solver.

## Default checks

The default termination check and convergence check functions are described below:
```
@doc
    Term_Check(x::BnBSolver,y::BnBModel,k_int::Int64)
    Conv_Check(x::BnBSolver,ubd::Float64,lbd::Float64)
```

Currently, the default is to never repeat a node.
```
@doc
Repeat_Node_Default(x::BnBSolver,y::BnBModel{Interval{T}}, Xin::Vector{Interval{T}},Xout::Vector{Interval{T}})
```

## Fathoming

By default, nodes are fathomed on value dominance.
```
@doc
    fathom!(y::BnBModel)
```

## Pre-processing and post-processing
By default, the pre-processing and post-processing functions simply return the
input Interval/MCInterval type vector and the prior feasibility value.
