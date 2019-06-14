"""
    default_storage!

Stores the two nodes to the stack.
"""
function default_storage!(x::Optimizer, y1::NodeBB, y2::NodeBB)
    x.maximum_node_id += 1; x.current_node_count += 1; x.stack[x.maximum_node_id] = y1
    x.maximum_node_id += 1; x.current_node_count += 1; x.stack[x.maximum_node_id] = y2
end

"""
    single_storage!

Stores the one nodes to the stack.
"""
function single_storage!(x::Optimizer, y::NodeBB)
    x.maximum_node_id += 1; x.current_node_count += 1; x.stack[x.maximum_node_id] = y
end
