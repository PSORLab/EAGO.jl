"""
    EAGO.fathom!(y::BnBModel)

Fathoms the stack by value dominance.
"""
function fathom!(y::BnBModel)
  loc = find(z->(z>y.UBDg),y.LBD)
  deleteat!(y.LBD,loc)
  deleteat!(y.UBD,loc)
  deleteat!(y.box,loc)
  deleteat!(y.id,loc)
end
