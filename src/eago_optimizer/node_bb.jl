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
    "Lower bound of problem solution on nodeBB"
    lower_bound::Float64
    "Upper bound of problem solution on nodeBB"
    upper_bound::Float64
    "Depth of node in B&B tree."
    depth::Int64
    "Unique id for each node."
    id::Int64
end

# Constructors
NodeBB() = NodeBB(Float64[], Float64[], -Inf, Inf, 0, 1)
NodeBB(x::NodeBB) = NodeBB(x.lower_variable_bounds, x.upper_variable_bounds,
                           x.lower_bound, x.upper_bound, x.depth, x.id)

# Copy utilities
Base.copy(x::NodeBB) = NodeBB(copy(x.lower_variable_bounds),
                              copy(x.upper_variable_bounds),
                              x.lower_bound, x.upper_bound, x.depth, x.id)

# using alternative name as to not interfere with ordering...
function uninitialized(x::NodeBB)
    flag = isempty(x.lower_variable_bounds)
    flag &= isempty(x.upper_variable_bounds)
    flag &= x.lower_bound === -Inf
    flag &= x.upper_bound === Inf
    flag &= x.depth === 0
    flag &= x.id === 1
    return flag
end

# Access functions for broadcasting data easily
lower_variable_bounds(x::NodeBB) = x.lower_variable_bounds
upper_variable_bounds(x::NodeBB) = x.upper_variable_bounds
lower_variable_bounds(x::NodeBB, id::Int64, nid::Int64) = x.lower_variable_bounds[id:nid]
upper_variable_bounds(x::NodeBB, id::Int64, nid::Int64) = x.upper_variable_bounds[id:nid]
lower_bound(x::NodeBB) = x.lower_bound
upper_bound(x::NodeBB) = x.upper_bound
depth(x::NodeBB) = x.depth

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
diam(x::NodeBB) = x.upper_variable_bounds - x.lower_variable_bounds
mid(x::NodeBB) = 0.5*(x.upper_variable_bounds + x.lower_variable_bounds)
