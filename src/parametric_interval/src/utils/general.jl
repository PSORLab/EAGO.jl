"""
    is_equal(X1::Vector{IntervalType},X2::Vector{IntervalType}, atol::Float64)

Returns true if X1 and X2 are equal to within tolerance atol in all dimensions.
"""
function is_equal(X1,X2::Vector{IntervalType}, atol::Float64)
  out::Bool = true
  for i=1:length(X1)
    if (abs(X1[i].lo-X2[i].lo)>=atol ||
        abs(X1[i].hi-X2[i].hi)>=atol )
        out = false
        break
    end
  end
  return out
end
