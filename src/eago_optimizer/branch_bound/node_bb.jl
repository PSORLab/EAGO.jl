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
mutable struct NodeBB
    lower_variable_bounds::Vector{Float64}
    upper_variable_bounds::Vector{Float64}
    lower_bound::Float64
    upper_bound::Float64
    depth::Int
    last_branch::Int
    branch_direction::Bool
end
NodeBB() = NodeBB(Float64[],Float64[],-Inf,Inf,0,-1,false)

# Access functions for broadcasting data easily
lower_variable_bounds(x::NodeBB) = x.lower_variable_bounds
upper_variable_bounds(x::NodeBB) = x.upper_variable_bounds
lower_variable_bounds(x::NodeBB,id::Int,nid::Int) = x.lower_variable_bounds[id:nid]
upper_variable_bounds(x::NodeBB,id::Int,nid::Int) = x.upper_variable_bounds[id:nid]
lower_bound(x::NodeBB) = x.lower_bound
upper_bound(x::NodeBB) = x.upper_bound
depth(x::NodeBB) = x.depth
last_branch(x::NodeBB) = x.last_branch
branch_direction(x::NodeBB) = x.branch_direction

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

"""
    NodeHistory

Stores historical information associated with solving the problem.
"""
mutable struct NodeHistory
    lower_count::Int                       # Number of Lower Bounding Problems Solved
    upper_count::Int                       # Number of Upper Bounding Problems Solved
    lower_bound::Dict{Int,Float64}         # Iteration --> Lower Bound
    upper_bound::Dict{Int,Float64}         # Iteration --> Upper Bound
    lower_time::Dict{Int,Float64}          # Iteration --> Lower Problem Run Time
    upper_time::Dict{Int,Float64}          # Iteration --> Upper Problem Run Time
    preprocess_time::Dict{Int,Float64}     # Iteration --> Preprocess Problem Run Time
    postprocess_time::Dict{Int,Float64}    # Iteration --> Postprocess Problem Run Time
    cut_count::Dict{Int,Int}               # Iteration --> Number of Cuts Used
    count::Dict{Int,Int}                  # Iteration --> Number of Nodes
end
NodeHistory() = NodeHistory(0,0,Dict{Int,Float64}(0 => -Inf),
                                Dict{Int,Float64}(0 => Inf),
                                Dict{Int,Float64}(0 => 0.0),
                                Dict{Int,Float64}(0 => 0.0),
                                Dict{Int,Float64}(0 => 0.0),
                                Dict{Int,Float64}(0 => 0.0),
                                Dict{Int,Int}(0 => 0),
                                Dict{Int,Int}(0 => 0))

NodeBB(x::NodeBB) = NodeBB(x.lower_variable_bounds, x.upper_variable_bounds,
              x.lower_bound, x.upper_bound, x.depth, x.last_branch,
              x.branch_direction)

function copy(x::NodeBB)
    return NodeBB(x.lower_variable_bounds, x.upper_variable_bounds,
                  x.lower_bound, x.upper_bound, x.depth, x.last_branch,
                  x.branch_direction)
end

get_lower_bound(x::NodeHistory) = [x.lower_bound[i-1] for i in 1:length(x.lower_bound)]
get_upper_bound(x::NodeHistory) = [x.upper_bound[i-1] for i in 1:length(x.upper_bound)]
get_lower_time(x::NodeHistory) = [x.lower_time[i-1] for i in 1:length(x.lower_time)]
get_upper_time(x::NodeHistory) = [x.upper_time[i-1] for i in 1:length(x.upper_time)]
get_preprocess_time(x::NodeHistory) = [x.preprocess_time[i-1] for i in 1:length(x.preprocess_time)]
get_postprocess_time(x::NodeHistory) = [x.postprocess_time[i-1] for i in 1:length(x.postprocess_time)]
get_solution_time(x::NodeHistory) = get_lower_time(x) + get_upper_time(x) + get_preprocess_time(x) + get_postprocess_time(x)
get_iteration_number(x::NodeHistory) = [i-1 for i in 1:length(x.count)]
get_node_count(x::NodeHistory) = [x.count[i-1] for i in 1:length(x.count)]
get_absolute_gap(x::NodeHistory) =  get_upper_bound(x) - get_lower_bound(x)
function get_relative_gap(x::NodeHistory)
    raw_relative_gap = abs.(get_upper_bound(x) - get_lower_bound(x))./max.(abs.(get_upper_bound(x)), abs.(get_lower_bound(x)))
    [isnan(raw_relative_gap[i]) ? Inf : raw_relative_gap[i] for i in 1:length(raw_relative_gap)]
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
