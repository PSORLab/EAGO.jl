"""
    extended_divide

Subfunction to generate output for extended division.
"""
function extended_divide(A::IntervalType)
  if ((A.lo == -0.0) && (A.hi == 0.0))
    B::IntervalType = Interval(-Inf,Inf)
    C::IntervalType = B
    return 0,B,C
  end
  if (A.lo == 0.0)
    B = Interval(1.0/A.hi,Inf)
    C = Interval(Inf,Inf)
    return 1,B,C
  elseif (A.hi == 0.0)
    B = Interval(-Inf,1.0/A.lo)
    C = Interval(-Inf,-Inf)
    return 2,B,C
  else
    B = Interval(-Inf,1.0/A.lo)
    C = Interval(1.0/A.hi,Inf)
    return 3,B,C
  end
end

"""
    extended_process

Generates output boxes for extended division and flag.
"""
# Error & Case Handling for Extended IA in NewtonGS (CHECKED C++)
function extended_process(N::IntervalType,X::IntervalType,Mii::IntervalType,
                          S1::IntervalType,S2::IntervalType,B::IntervalType,rtol::Float64)
  v = 1
  Ntemp::IntervalType = copy(N)
  M::IntervalType = (B+S1+S2)+Interval(-rtol,rtol)
  if (M.lo<=0 && M.hi>=0)
    return 0, Interval(-Inf,Inf), Ntemp
  end
  if (v == 1)
    k,IML::IntervalType,IMR::IntervalType = extended_divide(Mii)
    if (k == 1)
      return 0, (mid(X)-M*IML), Ntemp
    elseif (k == 2)
      return 0, (mid(X)-M*IMR), Ntemp
    elseif (k == 3)
      NR = mid(X)-M*IMR
      NL = mid(X)-M*IML
      if (~isdisjoint(NL,X) && isdisjoint(NR,X))
        return 0, NL, Ntemp
      elseif (~isdisjoint(NR,X) && isdisjoint(NL,X))
        return 0, NR, Ntemp
      elseif (~isdisjoint(NL,X) && ~isdisjoint(NR,X))
        N = NL
        Ntemp = NR
        return 1, NL, NR
      else
        return -1, N, Ntemp
      end
    end
  end
  return 0, N, Ntemp
end
