function sinh(a::MCInterval{T}) where T<:AbstractFloat
    isempty(a) && return a

    return MCInterval{T}(sinh(a.lo), sinh(a.hi))
end

function cosh(a::MCInterval{T}) where T<:AbstractFloat
    isempty(a) && return a

    return MCInterval{T}(cosh(mig(a)), cosh(mag(a)))
end

function tanh(a::MCInterval{T}) where T<:AbstractFloat
    isempty(a) && return a

    return MCInterval{T}(tanh(a.lo), tanh(a.hi))
end

function asinh(a::MCInterval{T}) where T<:AbstractFloat
    isempty(a) && return a

    return MCInterval{T}(asinh(a.lo), asinh(a.hi))
end

function acosh(a::MCInterval{T}) where T<:AbstractFloat
    domain = MCInterval(one(T), Inf)
    a = a ∩ domain
    isempty(a) && return a

    return MCInterval{T}(acosh(a.lo), acosh(a.hi))
end

function atanh(a::MCInterval{T}) where T<:AbstractFloat
    domain = MCInterval{T}(-one(T), one(T))
    a = a ∩ domain

    isempty(a) && return a

    res_lo = atanh(a.lo)
    res_hi = atanh(a.hi)

    (res_lo == res_hi == Inf || res_lo == res_hi == -Inf) && return emptyMCinterval(T)

    return MCInterval{T}(res_lo, res_hi)
end
