"""
    fathom!

Selects and deletes nodes from stack with lower bounds greater than global
upper bound.
"""
function fathom!(d::Optimizer)
  if ~isempty(d.stack)
    # Find the lowest upper bound without delecting
    MinKey, MinNode = first(d.stack)
    for (Key, Node) in d.stack
      if Node.upper_bound < d.global_upper_bound
        MinKey = Key
        MinNode = Node
        d.global_upper_bound = Node.upper_bound
      end
    end
    # Deletes all nodes with upper bound greater than minimum
    for (Key, Node) in d.stack
      if Node.lower_bound > d.global_upper_bound
          delete!(d.stack, Key)
      end
    end
  end
  d.current_node_count = length(d.stack)
end
