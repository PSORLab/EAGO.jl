function bisect(X::MCInterval)
    m = mid(X)
    return (MCInterval(X.lo, m), MCInterval(m, X.hi))
end

"""
    EAGO.Bisect_Abs(S::BnBSolver,B::BnBModel{T},N::Vector{T})

Returns two interval boxes 'X1,X2' created by bisecting 'N' in the highest width dimension.
"""
function Bisect_Abs(S::BnBSolver,B::BnBModel{T},N::Vector{T}) where {T}
  i::Int64 = indmax(diam.(N))
  N1::T,N2::T = bisect(N[i])
  X1::Vector{T} = deepcopy(N)
  X2::Vector{T} = deepcopy(N)
  X1[i] = N1
  X2[i] = N2
  return X1,X2
end

"""
    EAGO.Bisect_Rel(S::BnBSolver,B::BnBModel,N::Vector{Interval{Float64}})

Returns two interval boxes 'X1,X2' created by bisecting 'N' in the highest width dimension after scaling by initial box size.
"""
function Bisect_Rel(S::BnBSolver,B::BnBModel{T},N::Vector{T}) where {T}
  i::Int64 = indmax(diam.(N)./diam.(B.Init_Box))
  N1::T,N2::T = bisect(N[i])
  X1::Vector{T} = deepcopy(N)
  X2::Vector{T} = deepcopy(N)
  X1[i] = N1
  X2[i] = N2
  return X1,X2
end

"""
    EAGO.Bisect_Abs_Imp(S::BnBSolver,B::BnBModel,N::Vector{Interval{Float64}})

Returns two interval boxes 'X1,X2' created by bisecting 'N' in the highest width
dimension greater than 'nx'.
"""
function Bisect_Abs_Imp(S::BnBSolver,B::BnBModel{T},N::Vector{T},nx::Q) where {T,Q<:Integer}
  i::Int64 = indmax(diam.(N[(nx+1):end]))
  N1::T,N2::T = bisect(N[nx+i])
  X1::Vector{T} = deepcopy(N)
  X2::Vector{T} = deepcopy(N)
  X1[nx+i] = N1
  X2[nx+i] = N2
  return X1,X2
end

"""
    EAGO.Bisect_Rel_Imp(S::BnBSolver,B::BnBModel,N::Vector{Interval{Float64}})

Returns two interval boxes 'X1,X2' created by bisecting 'N' in the highest width
dimension greater than 'nx' after scaling by initial box size.
"""
function Bisect_Rel_Imp(S::BnBSolver,B::BnBModel{T},N::Vector{T},nx::Q) where {T,Q<:Integer}
  i::Int64 = indmax(diam.(N[(nx+1):end])./diam.(B.Init_Box[(nx+1):end]))
  N1::T,N2::T = bisect(N[nx+i])
  X1::Vector{T} = deepcopy(N)
  X2::Vector{T} = deepcopy(N)
  X1[nx+i] = N1
  X2[nx+i] = N2
  return X1,X2
end
