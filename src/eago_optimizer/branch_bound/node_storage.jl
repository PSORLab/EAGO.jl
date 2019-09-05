"""
    default_storage!

Stores the two nodes to the stack.
"""
function node_storage!(x::Optimizer, y1::NodeBB, y2::NodeBB)
    x.maximum_node_id += 1; x.current_node_count += 1; x.stack[x.maximum_node_id] = y1
    x.maximum_node_id += 1; x.current_node_count += 1; x.stack[x.maximum_node_id] = y2
end

"""
    single_storage!

Stores the one nodes to the stack.
"""
function single_storage!(x::Optimizer, y::NodeBB)
    x.maximum_node_id += 1
    x.current_node_count += 1
    lower_bound = max(y.lower_bound, x.current_lower_info.value)
    upper_bound = min(y.upper_bound, x.current_upper_info.value)
    n = NodeBB(y.lower_variable_bounds, y.upper_variable_bounds, lower_bound, upper_bound, y.depth + 1, -1, false)
    x.stack[x.maximum_node_id] = n
end
