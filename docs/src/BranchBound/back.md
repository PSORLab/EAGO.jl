## Default checks

Below are the default function for the Branch-and-Bound library. The EAGO solver
populates these fields based on user inputs to the solver in order to deliver a
valid nonconvex NLP solver.

The default termination check and convergence check functions are described below:
```
@docs
    EAGO.Term_Check(x::BnBSolver,y::BnBModel,k_int::Int64)
    EAGO.Conv_Check(x::BnBSolver,ubd::Float64,lbd::Float64)
```

Currently, the default is to never repeat a node.
```
@docs
    EAGO.Repeat_Node_Default(x::BnBSolver,y::BnBModel{Interval{T}}, Xin::Vector{Interval{T}},Xout::Vector{Interval{T}})
```

## Fathoming

By default, nodes are fathomed on value dominance.
```
@docs
    EAGO.fathom!(y::BnBModel)
```

## Pre-processing and post-processing
By default, the pre-processing and post-processing functions simply return the
input Interval/MCInterval type vector and the prior feasibility value.

## Bisection Methods
Method for absolute width bisection on all dimension in stack:

```@docs
    EAGO.Bisect_Abs
```

Method for relative width bisection on all dimension in stack:

```@docs
    EAGO.Bisect_Rel
```

Method for absolute width bisection ignore first `nx` dimensions:

```@docs
    EAGO.Bisect_Abs_Imp
```

Method for relative width bisection ignore first `nx` dimensions::

```@docs
    EAGO.Bisect_Rel_Imp
```
## Storage Methods for Common Branching Schemes
Node storage method for breadth-first search:

```@docs
    EAGO.BM_breadth!
```

Node storage method for depth-first or best-first search:

```@docs
    EAGO.BM_depth_best!
```

Node storage method for adding a single node to the top of the stack:

```@docs
    EAGO.BM_Single!
```

## Selection Methods for Common Branching Schemes

Select node for best-first search:

```@docs
    EAGO.NS_best(B::BnBModel)
```

Select node for depth-first or breadth-first search:

```@docs
    EAGO.NS_depth_breadth(B::BnBModel)
```

## Functions for generating console displayed
```@docs
    EAGO.print_int!(B::BnBSolver,k_int::Int64,k_nod::Int64,nid::Int64,lbdp::Float64,lbd::Float64, ubd::Float64,feasL::Bool,feasU::Bool)
```
```@docs
    EAGO.print_results!(B::BnBSolver,sol::Float64,pnt,feas::Bool,lbd_bool::Bool)
```
```@docs
    EAGO.print_sol!(x::BnBSolver,y::BnBModel,ubdcnt::Int64,lbdcnt::Int64,ubdtime::Float64,lbdtime::Float64)
```
