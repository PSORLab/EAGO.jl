"""
    EAGO.Strict_XinY(X::Vector{Interval{T}},Y::Vector{Interval{T}})

Returns true if X is strictly in Y (X.lo>Y.lo && X.hi<Y.hi).
"""
function Strict_XinY(X::Vector{Interval{T}},Y::Vector{Interval{T}}) where {T<:AbstractFloat}
  k::Bool = true
  for i=1:length(X)
    if ((X[i].lo<=Y[i].lo)||
        (X[i].hi>=Y[i].hi))
      k = false
    end
  end
  return k
end

"""
    EAGO.Strict_XinY(X::Vector{MCInterval{T}},Y::Vector{MCInterval{T}})

Returns true if X is strictly in Y (X.lo>Y.lo && X.hi<Y.hi).
"""
function Strict_XinY(X::Vector{MCInterval{T}},Y::Vector{MCInterval{T}}) where {T<:AbstractFloat}
  k::Bool = true
  for i=1:length(X)
    if ((X[i].lo<=Y[i].lo)||
        (X[i].hi>=Y[i].hi))
      k = false
    end
  end
  return k
end

"""
    Strict_XinY(X::Interval{T},Y::Interval{T})

Returns true if X is strictly in Y (X.lo>Y.lo && X.hi<Y.hi).
"""
function Strict_XinY(X::Interval{T},Y::Interval{T}) where {T<:AbstractFloat}
  (X.lo<=Y.lo) && return false
  (X.hi>=Y.hi) && return false
  return true
end

"""
    Strict_XinY(X::MCInterval{T},Y::MCInterval{T})

Returns true if X is strictly in Y (X.lo>Y.lo && X.hi<Y.hi).
"""
function Strict_XinY(X::MCInterval{T},Y::MCInterval{T}) where {T<:AbstractFloat}
  (X.lo<=Y.lo) && return false
  (X.hi>=Y.hi) && return false
  return true
end

"""
    Strictly_InRoutine

Checks that N is strictly in the interior of X. Sets the upper, lower, and both
inclusions storage functions as appropriate.
"""
function Strictly_InRoutine!(i::Integer,N::Vector{V},X::Vector{V},incl::Vector{Bool},
                            inclL::Vector{Bool},inclH::Vector{Bool}) where V<:Union{MCInterval,Interval}
    if Strict_XinY(N[i],X[i])
        incl[i] = true
        inclH[i] = true
        inclL[i] = true
    else
        incl[i] = false
        inclL[i] = false
        inclH[i] = false
        if (N[i].lo>X[i].lo)
            inclL[i] = true
        elseif (N[i].hi<X[i].hi)
            inclH[i] = true
        end
    end
end

"""
    InclTest(Iflg::Bool,inc::Vector{Bool},nx::Integer)

Takes the inclusion flag `Iflg` and sets it to true if all dimensions of the
state space were included `inc` is a vector of length `nx` containing a boolean
that indicates whether or not that dimension was included in the prior box. If
inclusions is satisfied for any box in the contractor it is true always.
"""
function InclTest(Iflg::Bool,inc::Vector{Bool},nx::Integer)
    if (Iflg == false)
      for i=1:nx
        if (inc[i] == true)
          Iflg = true
          continue
        elseif (inc[i] == false)
          Iflg = false
          break
        end
      end
    end
    return Iflg
end
