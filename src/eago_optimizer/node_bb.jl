"""
    NodeBB

Stores information associated with each node in Branch & Bound tree.
- `lower_variable_bounds::Vector{Float64}`: Lower bounds of variable box.
- `upper_variable_bounds::Vector{Float64}`: Upper bounds of variable box.
- `lower_bound::Float64`: Lower bound of problem solution on nodeBB
- `upper_bound::Float64`: Upper bound of problem solution on nodeBB
- `depth::Int`: Depth of node in B&B tree.
- `last_branch::Int`: Last dimension branched on.
- `branch_direction::Bool`: For future use with branching heurestics (false if nlp solve not in box, true otherwise).
"""
struct NodeBB
    lower_variable_bounds::Vector{Float64}
    upper_variable_bounds::Vector{Float64}
    lower_bound::Float64
    upper_bound::Float64
    depth::Int64
    id::Int64
end
NodeBB() = NodeBB(Float64[], Float64[], -Inf, Inf, 0, 1)
isless(x::NodeBB, y::NodeBB) = x.lower_bound < y.lower_bound

# Access functions for broadcasting data easily
lower_variable_bounds(x::NodeBB) = x.lower_variable_bounds
upper_variable_bounds(x::NodeBB) = x.upper_variable_bounds
lower_variable_bounds(x::NodeBB, id::Int64, nid::Int64) = x.lower_variable_bounds[id:nid]
upper_variable_bounds(x::NodeBB, id::Int64, nid::Int64) = x.upper_variable_bounds[id:nid]
lower_bound(x::NodeBB) = x.lower_bound
upper_bound(x::NodeBB) = x.upper_bound
depth(x::NodeBB) = x.depth
last_branch(x::NodeBB) = x.last_branch
branch_direction(x::NodeBB) = x.branch_direction

function diam(x::NodeBB)
    return x.upper_variable_bounds - x.lower_variable_bounds
end
function mid(x::NodeBB)
    return 0.5*(x.upper_variable_bounds + x.lower_variable_bounds)
end

length(x::NodeBB) = length(x.lower_variable_bounds)
function isempty(x::NodeBB)
    bool = false
    for i in 1:length(x)
        @inbounds lower = x.lower_variable_bounds[i]
        @inbounds upper = x.upper_variable_bounds[i]
        if (lower >= upper)
            bool = true
            break
        end
    end
    bool
end

NodeBB(x::NodeBB) = NodeBB(x.lower_variable_bounds, x.upper_variable_bounds,
                           x.lower_bound, x.upper_bound, x.depth, x.id)

function copy(x::NodeBB)
    return NodeBB(x.lower_variable_bounds, x.upper_variable_bounds,
                  x.lower_bound, x.upper_bound, x.depth, x.id)
end
"""
    same_box

Checks that node `x` and `y` have equal domains withing a tolerance of `atol`.
"""
function same_box(x::NodeBB,y::NodeBB, atol::Float64)
    bool = true
    (isempty(x.lower_variable_bounds) ⊻ isempty(y.lower_variable_bounds)) && (return false)
    (isempty(x.upper_variable_bounds) ⊻ isempty(y.upper_variable_bounds)) && (return false)
    for i in 1:length(x)
        if (x.lower_variable_bounds[i] != y.lower_variable_bounds[i])
            bool = false; break
        end
        if (x.upper_variable_bounds[i] != y.upper_variable_bounds[i])
            bool = false; break
        end
    end
    return bool
end
