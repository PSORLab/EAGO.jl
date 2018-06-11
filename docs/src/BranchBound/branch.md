## Storage Methods for Common Branching Schemes
Node storage method for breadth-first search:

```@docs
BM_breadth!(S::BnBSolver,B::BnBModel{T},tL::Float64,tU::Float64,X1::Vector{T},X2::Vector{T},pos::Int64)
```

Node storage method for depth-first or best-first search:

```@docs
BM_depth_best!(S::BnBSolver,B::BnBModel{T},tL::Float64,tU::Float64,X1::Vector{T}, X2::Vector{T},pos::Int64)
```

Node storage method for adding a single node to the top of the stack:

```@docs
BM_Single!(S::BnBSolver,B::BnBModel{T},tL::Float64,tU::Float64,X::Vector{T},pos::Int64)
```

## Selection Methods for Common Branching Schemes

Select node for best-first search:

```@docs
NS_best(B::BnBModel)
```

Select node for depth-first or breadth-first search:

```@docs
NS_depth_breadth(B::BnBModel)
```
