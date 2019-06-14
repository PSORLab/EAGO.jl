"""
    node_select_best!

Selects node with the lowest lower bound in stack.
"""
function node_select_best!(d::Optimizer)
  minkey, minnode = first(d.stack)
  minvalue = minnode.lower_bound
  for (key, node) in d.stack
    if node.lower_bound < minvalue
      minkey = key
      minnode = node
      minvalue = node.lower_bound
    end
  end
  delete!(d.stack,minkey)
  minkey,minnode
end
