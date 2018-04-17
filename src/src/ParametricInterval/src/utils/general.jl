"""
    isEqual(X1::Vector{Interval{T}},X2::Vector{Interval{T}},atol::Float64)

Returns true if X1 and X2 are equal to within tolerance atol in all dimensions.
"""
function isEqual(X1::Vector{Interval{T}},X2::Vector{Interval{T}},
                 atol::Float64) where {T<:AbstractFloat}
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
function isEqual(X1::Vector{MCInterval{T}},X2::Vector{MCInterval{T}},
                 atol::Float64) where {T<:AbstractFloat}
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
