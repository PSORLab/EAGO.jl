"""
    EAGO.NS_best(B::BnBModel)

Takes a single input `B::BnBModel`. Selects node with the lowest upper lower bound. Returns `(IntvBox,LBD,UBD,id,pos)`
where `Intv` is the interval box, `LBD` is the lower bound of the node, `UBD` is the
upper bound of the node, `id` is the id number of the node, and `pos` is the position
of the node in the BnB tree.
"""
function NS_best(B::BnBModel)
  tL,ind = findmin(B.LBD)
  return splice!(B.box,ind),splice!(B.LBD,ind),splice!(B.UBD,ind),splice!(B.id,ind),splice!(B.pos,ind)
end

"""
    EAGO.NS_depth_breadth(B::BnBModel)

Takes a single input `B::BnBModel`. Selects node on the top of the stack. Returns `(IntvBox,LBD,UBD,id,pos)`
where `Intv` is the intervalbox, `LBD` is the lower bound of the node, `UBD` is the
upper bound of the node, `id` is the id number of the node, and `pos` is the position
of the node in the BnB tree.
"""
function NS_depth_breadth(B::BnBModel)
  return pop!(B.box),pop!(B.LBD),pop!(B.UBD),pop!(B.id),pop!(B.pos)
end
