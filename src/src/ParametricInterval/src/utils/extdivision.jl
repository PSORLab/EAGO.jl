"""
    extDivide

Subfunction to generate output for extended division.
"""
function extDivide(A::Interval{T}) where {T<:AbstractFloat}
  #println("start extdivide")
  if ((A.lo == -0.0) && (A.hi == 0.0))
    B::Interval{T} = Interval(-Inf,Inf)
    C::Interval{T} = B
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
function extDivide(A::MCInterval{T}) where {T<:AbstractFloat}
  #println("start extdivide")
  if ((A.lo == -0.0) && (A.hi == 0.0))
    B::MCInterval{T} = MCInterval{T}(-Inf,Inf)
    C::MCInterval{T} = B
    return 0,B,C
  end
  if (A.lo == 0.0)
    B = MCInterval{T}(1.0/A.hi,Inf)
    C = MCInterval{T}(Inf,Inf)
    return 1,B,C
  elseif (A.hi == 0.0)
    B = MCInterval{T}(-Inf,1.0/A.lo)
    C = MCInterval{T}(-Inf,-Inf)
    return 2,B,C
  else
    B = MCInterval{T}(-Inf,1.0/A.lo)
    C = MCInterval{T}(1.0/A.hi,Inf)
    return 3,B,C
  end
end

"""
    extProcess

Generates output boxes for extended division and flag.
"""
# Error & Case Handling for Extended IA in NewtonGS (CHECKED C++)
function extProcess(N::Interval{T},X::Interval{T},Mii::Interval{T},
                    S1::Interval{T},S2::Interval{T},B::Interval{T},rtol::Float64) where {T<:AbstractFloat}
  v = 1
  Ntemp::Interval{T} = copy(N)
  M::Interval{T} = (B+S1+S2)+Interval(-rtol,rtol)
  if (M.lo<=0 && M.hi>=0)
    return 0, Interval(-Inf,Inf), Ntemp
  end
  if (v == 1)
    k,IML::Interval{T},IMR::Interval{T} = extDivide(Mii)
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

function extProcess(N::MCInterval{T},X::MCInterval{T},Mii::MCInterval{T},
                    S1::MCInterval{T},S2::MCInterval{T},B::MCInterval{T},rtol::Float64) where {T<:AbstractFloat}
  v = 1
  Ntemp::MCInterval{T} = copy(N)
  M::MCInterval{T} = (B+S1+S2)+MCInterval(-rtol,rtol)
  if (M.lo<=0 && M.hi>=0)
    return 0, MCInterval(-Inf,Inf), Ntemp
  end
  if (v == 1)
    k,IML::MCInterval{T},IMR::MCInterval{T} = extDivide(Mii)
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
