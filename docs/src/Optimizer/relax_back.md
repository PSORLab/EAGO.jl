# Nonlinear Backend

## Graphs, Caches, Forward and Reverse Propagation

EAGO makes use of a specialized tape structure for each function in order to compute valid
composite bounds and relaxations. Each variable, constant, and expression is respresented by
a node in a directed graph structure. 

```@docs
    EAGO.Node
    EAGO.NodeClass
    EAGO.AtomType
```

```@docs
    EAGO.DirectedTree
```

EAGO organizes information associated with each node a given graph structure using an
`EAGO.AbstractCache` which stores the given information.

```@docs
    EAGO.AbstractCache
```

Information in a given `EAGO.AbstractCache` is populated by performing a series of 
forward and reverse passes of the graph structure which dispatch off of an
`EAGO.AbstractCacheAttribute` which indicates what particular information is desired.

```@docs
    EAGO.AbstractCacheAttribute
```
Three included `AbstractCacheAttributes` are used to 

```@docs
    EAGO.Relax
    EAGO.RelaxAA
    EAGO.RelaxMulEnum
```

The forward and reverse routines are overloaded as follows:

XXX

## Other routines
```@docs
    EAGO.is_safe_cut!(m::GlobalOptimizer, f::MathOptInterface.ScalarAffineFunction{Float64}) 
```