# Copyright (c) 2018: Matthew Wilhelm & Matthew Stuber.
# This code is licensed under MIT license (see LICENSE.md for full details)
#############################################################################
# EAGO
# A development environment for robust and global optimization
# See https://github.com/PSORLab/EAGO.jl
#############################################################################
# src/eago_optimizer/node_bb.jl
# Defines storage for a node in the B&B tree & utilities functions
#############################################################################

@enum(BranchDirection, BD_NONE, BD_NEG, BD_POS)

# Used internally to set & get variables in full problem space or just branch variables
struct FullVar end
struct BranchVar end

"""
$(TYPEDEF)

Stores information associated with each node in Branch & Bound tree.

$(TYPEDFIELDS)
"""
struct NodeBB
    "Lower bounds of variable box."
    lower_variable_bounds::Vector{Float64}
    "Upper bounds of variable box."
    upper_variable_bounds::Vector{Float64}
    "Is dimension integer valued"
    is_integer::BitVector
    "Are all dimensions continuous or fixed"
    continuous::Bool
    "Lower bound of problem solution on nodeBB"
    lower_bound::Float64
    "Upper bound of problem solution on nodeBB"
    upper_bound::Float64
    "Depth of node in B&B tree."
    depth::Int
    "Unique id for each node."
    id::Int
    "Whether last branch was negative or positive in direction"
    branch_direction::BranchDirection
    "Dimension of last branch"
    last_branch::Int
    "Extent of last branch (using for psuedocost calculation)"
    branch_extent::Float64
end

# Constructors
function NodeBB(l::Vector{Float64}, u::Vector{Float64}, d::BitVector)
    NodeBB(l, u, d, any(d), -Inf, Inf, 1, 1, BD_NONE, -1, 0.0)
end
NodeBB() = NodeBB(Float64[], Float64[], BitVector(), false, -Inf, Inf, 0, 1, BD_NONE, -1, 0.0)
NodeBB(x::NodeBB) = NodeBB(copy(x.lower_variable_bounds), copy(x.upper_variable_bounds),
                           copy(x.is_integer), x.continuous,
                           x.lower_bound, x.upper_bound, x.depth, x.id,
                           x.branch_direction, x.last_branch, x.branch_extent)

# Copy utilities
Base.copy(x::NodeBB) = NodeBB(copy(x.lower_variable_bounds),
                              copy(x.upper_variable_bounds),
                              copy(x.is_integer), x.continuous,
                              x.lower_bound, x.upper_bound, x.depth, x.id,
                              x.branch_direction, x.last_branch, x.branch_extent)

# using alternative name as to not interfere with ordering...
function uninitialized(x::NodeBB)
    flag = isempty(x.lower_variable_bounds)
    flag &= isempty(x.upper_variable_bounds)
    flag &= isempty(x.is_integer)
    flag &= x.lower_bound === -Inf
    flag &= x.upper_bound === Inf
    flag &= x.depth === 0
    flag &= x.id === 1
    return flag
end

# Access functions for broadcasting data easily
@inline lower_variable_bounds(x::NodeBB) = x.lower_variable_bounds
@inline upper_variable_bounds(x::NodeBB) = x.upper_variable_bounds
@inline lower_variable_bounds(x::NodeBB, i::Int) = x.lower_variable_bounds[i]
@inline upper_variable_bounds(x::NodeBB, i::Int) = x.upper_variable_bounds[i]
@inline lower_variable_bounds(x::NodeBB, id::Int, nid::Int) = x.lower_variable_bounds[id:nid]
@inline upper_variable_bounds(x::NodeBB, id::Int, nid::Int) = x.upper_variable_bounds[id:nid]
@inline is_integer(x::NodeBB) = x.is_integer
@inline is_integer(x::NodeBB, id::Int) = x.is_integer[id]
@inline lower_bound(x::NodeBB) = x.lower_bound
@inline upper_bound(x::NodeBB) = x.upper_bound
@inline depth(x::NodeBB) = x.depth

# Iterations Functions
Base.isless(x::NodeBB, y::NodeBB) = x.lower_bound < y.lower_bound
Base.length(x::NodeBB) = length(x.lower_variable_bounds)
function Base.isempty(x::NodeBB)
    for i = 1:length(x)
        @inbounds lower = x.lower_variable_bounds[i]
        @inbounds upper = x.upper_variable_bounds[i]
        (lower > upper) && (return true)
    end
    return false
end

"""
$(FUNCTIONNAME)

Checks that node `x` and `y` have equal domains withing a tolerance of `atol`.
"""
function same_box(x::NodeBB, y::NodeBB, r::Float64)
    (isempty(x.lower_variable_bounds) ⊻ isempty(y.lower_variable_bounds)) && (return false)
    (isempty(x.upper_variable_bounds) ⊻ isempty(y.upper_variable_bounds)) && (return false)
    for i = 1:length(x)
        ~isapprox(x.lower_variable_bounds[i], y.lower_variable_bounds[i], atol=r) && (return false)
        ~isapprox(x.upper_variable_bounds[i], y.upper_variable_bounds[i], atol=r) && (return false)
    end
    return true
end

# Compute middle & diameter
@inline diam(x::NodeBB) = upper_variable_bounds(x) - lower_variable_bounds(x)
@inline diam(x::NodeBB, i::Int) = upper_variable_bounds(x,i) - lower_variable_bounds(x,i)

@inline mid(x::NodeBB) = 0.5*(upper_variable_bounds(x) + lower_variable_bounds(x))
@inline mid(x::NodeBB, i::Int) = 0.5*(upper_variable_bounds(x,i) + lower_variable_bounds(x,i))
