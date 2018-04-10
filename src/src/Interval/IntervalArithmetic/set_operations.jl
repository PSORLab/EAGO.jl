function in(x::T, a::MCInterval{T}) where T<:AbstractFloat
    isinf(x) && return false
    a.lo <= x <= a.hi
end

function isinterior(a::MCInterval{T}, b::MCInterval{T}) where T<:AbstractFloat
    isempty(a) && return true
    islessprime(b.lo, a.lo) && islessprime(a.hi, b.hi)
end

function isdisjoint(a::MCInterval{T}, b::MCInterval{T}) where T<:AbstractFloat
    (isempty(a) || isempty(b)) && return true
    islessprime(b.hi, a.lo) || islessprime(a.hi, b.lo)
end

function intersect(a::MCInterval{T}, b::MCInterval{T}) where T<:AbstractFloat
    isdisjoint(a,b) && return emptyMCinterval(T)

    MCInterval{T}(max(a.lo, b.lo), min(a.hi, b.hi))
end

intersect(a::MCInterval{T}, b::MCInterval{S}) where {T<:AbstractFloat,S<:AbstractFloat} = intersect(promote(a, b)...)
hull(a::MCInterval{T}, b::MCInterval{T}) where T<:AbstractFloat = MCInterval(min(a.lo, b.lo), max(a.hi, b.hi))
union(a::MCInterval{T}, b::MCInterval{T}) where T<:AbstractFloat = hull(a, b)

function setdiff(x::MCInterval{T}, y::MCInterval{T}) where T<:AbstractFloat
    intersection = x ∩ y

    isempty(intersection) && return [x]
    intersection == x && return typeof(x)[]  # x is subset of y; setdiff is empty

    x.lo == intersection.lo && return [MCInterval{T}(intersection.hi, x.hi)]
    x.hi == intersection.hi && return [MCInterval{T}(x.lo, intersection.lo)]

    return [MCInterval{T}(x.lo, y.lo),
            MCInterval{T}(y.hi, x.hi)]

end
