"""
    find_lower_bound

Selects node with the lowest lower bound in stack.
"""
function find_lower_bound(d::Optimizer)
  minkey, minnode = first(d.stack)
  minvalue = minnode.lower_bound
  for (key, node) in d.stack
    if node.lower_bound < minvalue
      minkey = key
      minnode = node
      minvalue = node.lower_bound
    end
  end
  minvalue
end
