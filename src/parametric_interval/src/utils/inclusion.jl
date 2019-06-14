"""
    strict_x_in_y(X::Vector{IntervalType},Y::Vector{IntervalType})

Returns true if X is strictly in Y (X.lo>Y.lo && X.hi<Y.hi).
"""
function strict_x_in_y(X::Vector{IntervalType},Y::Vector{IntervalType})
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
    strict_x_in_y(X::IntervalType, Y::IntervalType)

Returns true if X is strictly in Y (X.lo>Y.lo && X.hi<Y.hi).
"""
function strict_x_in_y(X::IntervalType, Y::IntervalType)
  (X.lo<=Y.lo) && return false
  (X.hi>=Y.hi) && return false
  return true
end

"""
    strictly_in_routine!

Checks that N is strictly in the interior of X. Sets the upper, lower, and both
inclusions storage functions as appropriate.
"""
function strictly_in_routine!(i::Int, N::Vector{IntervalType}, X,
                             incl::Vector{Bool}, inclL::Vector{Bool}, inclH::Vector{Bool})
    if strict_x_in_y(N[i],X[i])
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
    inclusion_test(inclusion_flag::Bool, inclusion_vector::Vector{Bool}, nx::Int)

Takes the inclusion flag `inclusion_flag` and sets it to true if all dimensions of the
state space were included `inclusion_vector` is a vector of length `nx` containing a boolean
that indicates whether or not that dimension was included in the prior box. If
inclusions is satisfied for any box in the contractor it is true always.
"""
function inclusion_test(inclusion_flag::Bool, inclusion_vector::Vector{Bool}, nx::Int)
    if (inclusion_flag == false)
      for i=1:nx
        if (inclusion_vector[i] == true)
          inclusion_flag = true
          continue
        elseif (inclusion_vector[i] == false)
          inclusion_flag = false
          break
        end
      end
    end
    return inclusion_flag
end
