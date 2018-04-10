# This file is part of the IntervalArithmetic.jl package; MIT licensed

function infty(T)
    (T == Float64) && return Inf64
    (T == Float32) && return Inf32
    return Inf16
end

function ninfty(T)
    (T == Float64) && return -Inf64
    (T == Float32) && return -Inf32
    return -Inf16
end

## Comparisons
function ==(a::MCInterval, b::MCInterval)
    isempty(a) && isempty(b) && return true
    a.lo == b.lo && a.hi == b.hi
end
!=(a::MCInterval, b::MCInterval) = !(a==b)


# Auxiliary functions: equivalent to </<=, but Inf <,<= Inf returning true
function islessprime(a::T, b::T) where T<:AbstractFloat
    (isinf(a) || isinf(b)) && a==b && return true
    a < b
end

# Weakly less, \le, <=
function <=(a::MCInterval, b::MCInterval)
    isempty(a) && isempty(b) && return true
    (isempty(a) || isempty(b)) && return false
    (a.lo ≤ b.lo) && (a.hi ≤ b.hi)
end

# Strict less: <
function <(a::MCInterval, b::MCInterval)
    isempty(a) && isempty(b) && return true
    (isempty(a) || isempty(b)) && return false
    islessprime(a.lo, b.lo) && islessprime(a.hi, b.hi)
end

# precedes
function precedes(a::MCInterval, b::MCInterval)
    (isempty(a) || isempty(b)) && return true
    a.hi ≤ b.lo
end

# strictpreceds
function strictprecedes(a::MCInterval, b::MCInterval)
    (isempty(a) || isempty(b)) && return true
    # islessprime(a.hi, b.lo)
    a.hi < b.lo
end
const ≺ = strictprecedes # \prec


# zero, one
zero(a::MCInterval{T}) where T<:AbstractFloat = MCInterval(zero(T))
zero(::Type{MCInterval{T}}) where T<:AbstractFloat = MCInterval(zero(T))
one(a::MCInterval{T}) where T<:AbstractFloat = MCInterval(one(T))
one(::Type{MCInterval{T}}) where T<:AbstractFloat = MCInterval(one(T))


## Addition and subtraction

+(a::MCInterval) = a
-(a::MCInterval) = MCInterval(-a.hi, -a.lo)

function +(a::MCInterval{T}, b::MCInterval{T}) where T<:AbstractFloat
    (isempty(a) || isempty(b)) && return emptyMCinterval(T)
    MCInterval{T}(a.lo + b.lo, a.hi + b.hi)
end

function -(a::MCInterval{T}, b::MCInterval{T}) where T<:AbstractFloat
    (isempty(a) || isempty(b)) && return emptyMCinterval(T)
    MCInterval{T}(a.lo - b.hi, a.hi - b.lo)
end


## Multiplication

function *(a::MCInterval{T}, b::MCInterval{T}) where T<:AbstractFloat
    (isempty(a) || isempty(b)) && return emptyMCinterval(T)

    (iszero(a) || iszero(b)) && return zero(MCInterval{T})

    if b.lo >= zero(T)
        a.lo >= zero(T) && return MCInterval{T}(a.lo*b.lo, a.hi*b.hi)
        a.hi <= zero(T) && return MCInterval{T}(a.lo*b.hi, a.hi*b.lo)
        return MCInterval{T}(a.lo*b.hi, a.hi*b.hi)   # zero(T) ∈ a
    elseif b.hi <= zero(T)
        a.lo >= zero(T) && return MCInterval{T}(a.hi*b.lo, a.lo*b.hi)
        a.hi <= zero(T) && return MCInterval{T}(a.hi*b.hi, a.lo*b.lo)
        return MCInterval{T}(a.hi*b.lo, a.lo*b.lo)   # zero(T) ∈ a
    else
        a.lo > zero(T) && return MCInterval{T}(a.hi*b.lo, a.hi*b.hi)
        a.hi < zero(T) && return MCInterval{T}(a.lo*b.hi, a.lo*b.lo)
        return MCInterval{T}(min(a.lo*b.hi, a.hi*b.lo), max(a.lo*b.lo, a.hi*b.hi))
    end
end


## Division

function inv(a::MCInterval{T}) where T<:AbstractFloat
    isempty(a) && return emptyMCinterval(T)

    if zero(T) ∈ a
        a.lo < zero(T) == a.hi && return MCInterval{T}(-Inf, inv(a.lo))
        a.lo == zero(T) < a.hi && return MCInterval{T}((inv(a.hi), Inf))
        a.lo < zero(T) < a.hi && return entireMCinterval(T)
        a == zero(a) && return emptyMCinterval(T)
    end

    MCInterval{T}(inv(a.hi), inv(a.lo))
end

function /(a::MCInterval{T}, b::MCInterval{T}) where T<:AbstractFloat

    S = typeof(a.lo / b.lo)
    (isempty(a) || isempty(b)) && return emptyMCinterval(S)
    iszero(b) && return emptyMCinterval(S)

    if b.lo > zero(T) # b strictly positive

        a.lo >= zero(T) && return MCInterval{T}(a.lo/b.hi, a.hi/b.lo)
        a.hi <= zero(T) && return MCInterval{T}(a.lo/b.lo, a.hi/b.hi)
        return MCInterval{T}(a.lo/b.lo, a.hi/b.lo)  # zero(T) ∈ a

    elseif b.hi < zero(T) # b strictly negative

        a.lo >= zero(T) && return MCInterval{T}(a.hi/b.hi, a.lo/b.lo)
        a.hi <= zero(T) && return MCInterval{T}(a.hi/b.lo, a.lo/b.hi)
        return MCInterval{T}(a.hi/b.hi, a.lo/b.hi)  # zero(T) ∈ a

    else   # b contains zero, but is not zero(b)

        iszero(a) && return zero(MCInterval{S})

        if iszero(b.lo)

            a.lo >= zero(T) && return MCInterval{T}(a.lo/b.hi, Inf)
            a.hi <= zero(T) && return MCInterval{T}(-Inf, a.hi/b.hi)
            return entireMCinterval(S)

        elseif iszero(b.hi)

            a.lo >= zero(T) && return MCInterval{T}(-Inf, a.lo/b.lo)
            a.hi <= zero(T) && return MCInterval{T}(a.hi/b.lo, Inf)
            return entireMCinterval(S)

        else

            return entireMCinterval(S)

        end
    end
end

//(a::MCInterval, b::MCInterval) = a / b    # to deal with rationals


function min_ignore_nans(args...)
    min(Iterators.filter(x->!isnan(x), args)...)
end

function max_ignore_nans(args...)
    max(Iterators.filter(x->!isnan(x), args)...)
end

## fma: fused multiply-add
function fma(a::MCInterval{T}, b::MCInterval{T}, c::MCInterval{T}) where T
    #T = promote_type(eltype(a), eltype(b), eltype(c))

    (isempty(a) || isempty(b) || isempty(c)) && return emptyMCinterval(T)

    if isentire(a)
        b == zero(b) && return c
        return entireMCinterval(T)

    elseif isentire(b)
        a == zero(a) && return c
        return entireMCinterval(T)

    end

    lo = min_ignore_nans(fma(a.lo, b.lo, c.lo),
                         fma(a.lo, b.hi, c.lo),
                         fma(a.hi, b.lo, c.lo),
                         fma(a.hi, b.hi, c.lo))

    hi = max_ignore_nans(fma(a.lo, b.lo, c.hi),
                         fma(a.lo, b.hi, c.hi),
                         fma(a.hi, b.lo, c.hi),
                         fma(a.hi, b.hi, c.hi))

    MCInterval{T}(lo, hi)
end


## Scalar functions on intervals (no directed rounding used)

function mag(a::MCInterval{T}) where T<:AbstractFloat
    isempty(a) && return convert(eltype(a), NaN)
    # r1, r2 = setrounding(T, RoundUp) do
    #     abs(a.lo), abs(a.hi)
    # end
    max( abs(a.lo), abs(a.hi) )
end

function mig(a::MCInterval{T}) where T<:AbstractFloat
    isempty(a) && return convert(eltype(a), NaN)
    zero(a.lo) ∈ a && return zero(a.lo)
    min(abs(a.lo), abs(a.hi))
end


# Infimum and supremum of an interval
inf(a::MCInterval) = a.lo
sup(a::MCInterval) = a.hi


## Functions needed for generic linear algebra routines to work
real(a::MCInterval) = a

function abs(a::MCInterval{T}) where {T<:AbstractFloat}
    isempty(a) && return emptyMCinterval(T)
    MCInterval{T}(mig(a), mag(a))
end

function min(a::MCInterval{T}, b::MCInterval{T}) where {T<:AbstractFloat}
    (isempty(a) || isempty(b)) && return emptyMCinterval(T)
    MCInterval{T}( min(a.lo, b.lo), min(a.hi, b.hi))
end

function max(a::MCInterval{T}, b::MCInterval{T}) where {T<:AbstractFloat}
    (isempty(a) || isempty(b)) && return emptyMCinterval(T)
    MCInterval{T}( max(a.lo, b.lo), max(a.hi, b.hi))
end

dist(a::MCInterval, b::MCInterval) = max(abs(a.lo-b.lo), abs(a.hi-b.hi))
eps(a::MCInterval) = max(eps(a.lo), eps(a.hi))

## floor, ceil, trunc, sign, roundTiesToEven, roundTiesToAway
function floor(a::MCInterval{T}) where {T<:AbstractFloat}
    isempty(a) && return emptyMCinterval(a)
    MCInterval{T}(floor(a.lo), floor(a.hi))
end

function ceil(a::MCInterval{T}) where {T<:AbstractFloat}
    isempty(a) && return emptyMCinterval(T)
    MCInterval{T}(ceil(a.lo), ceil(a.hi))
end

function trunc(a::MCInterval{T}) where {T<:AbstractFloat}
    isempty(a) && return emptyMCinterval(T)
    MCInterval{T}(trunc(a.lo), trunc(a.hi))
end

function sign(a::MCInterval{T}) where {T<:AbstractFloat}
    isempty(a) && return emptyMCinterval(T)
    return MCInterval{T}(sign(a.lo), sign(a.hi))
end

# mid, diam, radius
function mid(a::MCInterval{T},α) where {T<:AbstractFloat}

    isempty(a) && return convert(T, NaN)
    isentire(a) && return zero(a.lo)

    a.lo == -∞ && return nextfloat(-∞)
    a.hi == +∞ && return prevfloat(+∞)

    # @assert 0 ≤ α ≤ 1

    # return (1-α) * a.lo + α * a.hi  # rounds to nearest
    return α*(a.hi - a.lo) + a.lo  # rounds to nearest
end

function mid(a::MCInterval{T}) where {T<:AbstractFloat}

    isempty(a) && return convert(T, NaN)
    isentire(a) && return zero(a.lo)

    a.lo == -∞ && return nextfloat(a.lo)
    a.hi == +∞ && return prevfloat(a.hi)

    # @assert 0 ≤ α ≤ 1

    return 0.5 * (a.lo + a.hi)
end

function diam(a::MCInterval{T}) where {T<:AbstractFloat}
    isempty(a) && return convert(T, NaN)
    a.hi - a.lo
end

# Should `radius` this yield diam(a)/2? This affects other functions!
function radius(a::MCInterval{T}) where {T<:AbstractFloat}
    isempty(a) && return convert(eltype(a), NaN)
    m = mid(a)
    return max(m - a.lo, a.hi - m)
end

function step(x::MCInterval{T}) where {T}
      isempty(x) && return emptyMCinterval(T)
      xmin::T = ((x.lo)<zero(T)) ? zero(T) : one(T)
      xmax::T = ((x.hi)>=zero(T)) ? one(T) : zero(T)
      return MCInterval{T}(xmin,xmax)
end
