# Nonlinear Backend

## Graphs, Caches, Forward and Reverse Propagation

EAGO makes use of a specialized tape structure for each function in order to compute valid composite bounds and relaxations. Each variable, constant, and expression is respresented by a node in a directed graph structure. 

```@docs
    EAGO.Node
    EAGO.NodeClass
    EAGO.AtomType
```

```@docs
    EAGO.AbstractDirectedGraph
    EAGO.DirectedTree
```

Each field of the ith `EAGO.Node` using a basic access function. For instance the ith node's `ex_type` for graph `d` may be accessed by `ex_type(d, i)`. The `sparsity` of node `i` returns an ordered list of `children` nodes which form the argument tuple for the operator performed at node `i`. The `sparsity` of node `i` returns a list of `parent` nodes which form the argument tuple for the operator performed at node `i`. The `parameter_values` and `constant_values` functions are used to access the ith parameter values or ith constant values.

EAGO organizes information associated with each node in a given graph structure using an `EAGO.AbstractCache` which stores the given information.

```@docs
    EAGO.AbstractCache
    EAGO.initialize!(::AbstractCache, ::AbstractDirectedGraph)
```

Information in a given `EAGO.AbstractCache` is populated by performing a series of forward and reverse passes of the graph structure which dispatch off of an
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

```@docs
    EAGO.f_init!(t::AbstractCacheAttribute, g::AbstractDirectedGraph, c::AbstractCache)
    EAGO.fprop!(t::AbstractCacheAttribute, v::Variable, g::AbstractDirectedGraph, c::AbstractCache, k::Int)
    EAGO.fprop!(t::AbstractCacheAttribute, v::Subexpression, g::AbstractDirectedGraph, c::AbstractCache, k::Int)
    EAGO.fprop!(t::AbstractCacheAttribute, v::Expression, g::AbstractDirectedGraph, c::AbstractCache, k::Int)
    EAGO.fprop!(t::AbstractCacheAttribute, v::Parameter, g::AbstractDirectedGraph, c::AbstractCache, k::Int)
    EAGO.fprop!(t::AbstractCacheAttribute, v::Constant, g::AbstractDirectedGraph, c::AbstractCache, k::Int)
```

```@docs
    EAGO.r_init!(t::AbstractCacheAttribute, g::AbstractDirectedGraph, c::AbstractCache)
    EAGO.rprop!(t::AbstractCacheAttribute, v::Variable, g::AbstractDirectedGraph, c::AbstractCache, k::Int)
    EAGO.rprop!(t::AbstractCacheAttribute, v::Subexpression, g::AbstractDirectedGraph, c::AbstractCache, k::Int)
    EAGO.rprop!(t::AbstractCacheAttribute, v::Expression, g::AbstractDirectedGraph, c::AbstractCache, k::Int)
    EAGO.rprop!(t::AbstractCacheAttribute, v::Parameter, g::AbstractDirectedGraph, c::AbstractCache, k::Int)
    EAGO.rprop!(t::AbstractCacheAttribute, v::Constant, g::AbstractDirectedGraph, c::AbstractCache, k::Int)
```

Forward and reverse subroutines are overloaded for individual operators using through functions of the forms 
`fprop!(t::AbstractCacheAttribute, v::Val{AtomType}, g::AbstractDirectedGraph, b::AbstractCache, k::Int)` and
`rprop!(t::AbstractCacheAttribute, v::Val{AtomType}, g::AbstractDirectedGraph, b::AbstractCache, k::Int)`. 

## Other routines
```@docs
    EAGO.is_safe_cut!(m::GlobalOptimizer, f::MathOptInterface.ScalarAffineFunction{Float64}) 
```