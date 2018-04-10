"""
    tol_eq(X1,X2,atol)

Checks that values `X1` and `X2` are equal to within tolerance `atol`.
"""
function tol_eq(X1,X2,atol)
  out::Bool = true
  if (abs(X1-X2)>atol)
    out = false
  end
  return out
end
