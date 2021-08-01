
const USE_MIN_RANGE = true

struct AffineEAGO{N}
    c::Float64               # mid-point
    γ::SVector{N,Float64}    # affine terms
    Δ::Float64               # error term
end

function AffineEAGO(x::AffineEAGO{N}, p::Float64, q::Float64, δ::Float64) where N
    c = p * x.c + q
    γ = p .* x.γ
    Δ = p * x.Δ + δ
    AffineEAGO{N}(c, γ, δ)
end
AffineEAGO{N}(x::Float64, X::Interval{Float64}, i::Int) where N = AffineEAGO{N}(mid(X), seed_gradient(i, Val(N)), rad(X))

const UNIT_INTERVAL = Interval{Float64}(-1,1)
Interval(x::AffineEAGO{N}) where N = x.c + x.Δ*UNIT_INTERVAL + sum(y -> abs(y)*UNIT_INTERVAL, x.γ)
function bounds(x::AffineEAGO{N}) where N
    z = interval(x)
    z.lo, z.hi
end

zero(::Type{AffineEAGO{N}}) where N = AffineEAGO{N}(0.0, zero(SVector{N,Float64}), 0.0)
zero(x::AffineEAGO{N}) where N = AffineEAGO{N}(0.0, zero(SVector{N,Float64}), 0.0)

+(x::AffineEAGO{N}, y::AffineEAGO{N}) where N = AffineEAGO{N}(x.c + y.c, x.γ .+ y.γ, x.Δ + y.Δ)
+(x::AffineEAGO{N}, α::Float64) where N = AffineEAGO{N}(α+x.c, x.γ, x.Δ)
+(α::Float64, x::AffineEAGO{N}) where N = x + α

function *(x::AffineEAGO{N}, y::AffineEAGO{N}) where N
    c = x.c * y.c
    γ = x.c .* y.γ + y.c .* x.γ
    pos_sum = 0.0
    neg_sum = 0.0
    for i = 1:N
        vt = γ[i]*δ[i]
        if vt > 0.0
            pos_sum += vt
        else
            neg_sum -= vt
        end 
    end
    Δ = max(pos_sum, neg_sum)
    for i = 1:N
        for j = 1:(i-1)
            Δ += γ[i]*δ[j] + γ[j]*δ[i]
        end
    end
    return AffineEAGO{N}(c, γ, Δ)
end
*(x::AffineEAGO{N}, α::Float64) where N = AffineEAGO{N}(α*x.c, α.*x.γ, abs(α)*x.Δ)
*(α::Float64, x::AffineEAGO{N}) where N = x*α

function -(x::AffineEAGO{N}, y::AffineEAGO{N}) where N
    AffineEAGO{N}(x.c - y.c, x.γ .- y.γ, x.Δ + y.Δ)
end
-(x::AffineEAGO{N}) where N = AffineEAGO{N}(-x.c, .-(x.γ), x.Δ)
-(x::AffineEAGO{N}, α::Float64) where N = AffineEAGO{N}(x.c - α, x.γ, x.Δ)
-(α::Float64, x::AffineEAGO{N}) where N = α + (-x)

/(x::AffineEAGO{N}, α::Float64) where N = AffineEAGO{N}(x.c/α, x.γ/α, x.Δ/abs(α))
/(α::Float64, x::AffineEAGO{N}) where N = α*inv(x)

function log(x::AffineEAGO{N}) where N
    a, b = bounds(x)
    fa = log(a)
    fb = log(b)
    if USE_MIN_RANGE
        p = 1/b
        q = 0.5*(fa + fb - p*(a + b))
        Δ = abs(0.5*(fb - fa - p*(b - a)))
        return AffineEAGO(x, p, q, Δ)
    end
    p = (fb - fa)/(b - a)
    ξ = 1/p
    fξ = log(ξ)
    q = 0.5*(fa + fξ - p*(a + ξ))
    Δ = abs(0.5*(fξ - fa - p*(ξ - a)))
    return AffineEAGO(x, p, q, Δ)
end

function log10(x::AffineEAGO{N}) where N
    a, b = bounds(x)
    fa = log10(a)
    fb = log10(b)
    if USE_MIN_RANGE
        p = 1/(b*log(10))
        q = 0.5*(fa + fb - p*(a + b))
        Δ = abs(0.5*(fb - fa - p*(b - a)))
        return AffineEAGO(x, p, q, Δ)
    end
    p = (fb - fa)/(b - a)
    ξ = 1/p
    fξ = log10(ξ)
    q = 0.5*(fa + fξ - p*(a + ξ))
    Δ = abs(0.5*(fξ - fa - p*(ξ - a)))
    return AffineEAGO(x, p, q, Δ)
end

function pow_1d(x::AffineEAGO{N}, n::Int, p) where N
    a, b = bounds(x)
    fa = a^n
    fb = b^n
    if USE_MIN_RANGE
        q = 0.5*(fa + fb - p*(a + b))
        Δ = abs(0.5*(fb - fa - p*(b - a)))
        return AffineEAGO(x, p, q, Δ)
    end
    p = (fb - fa)/(b - a)
    ξ = (p/n)^(1/(n - 1))
    fξ = ξ^n
    q = 0.5*(fa + fξ - p*(a + ξ))
    Δ = abs(0.5*(fξ - fa - p*(ξ - a)))
    return AffineEAGO(x, p, q, Δ)
end

function pow_even(x::AffineEAGO{N}, n::Int) where N
    a, b = bounds(x)
    fa = a^n
    fb = b^n
    if USE_MIN_RANGE
        m = min(0.0, a, b)
        M = max(0.0, a, b)
        p = 0.0
        q = 0.5*(m + M)
        Δ = 0.5*(M - m)
        return AffineEAGO(x, p, q, Δ)
    end
    p = (fb - fa)/(b - a)
    ξ = (p/n)^(1/(n - 1))
    fξ = ξ^n
    q = 0.5*(fa + fξ - p*(a + ξ))
    Δ = abs(0.5*(fξ - fa - p*(ξ - a)))
    return AffineEAGO(x, p, q, Δ)
end

function pow_odd(x::AffineEAGO{N}, n::Int) where N
    a, b = bounds(x)
    fa = a^n
    fb = b^n
    p = (fb - fa)/(b - a)
    q = 0.0
    ξ = (p/n)^(1/(n-1))
    fξ = ξ^n
    Δ = abs(fξ - p*ξ)
    return AffineEAGO(x, p, q, Δ)
end

function ^(x::AffineEAGO{N}, n::Int) where N

    iszero(n) && zero(x)
    isone(n)  && one(x)

    xL, xU = bounds(x)
    if (xL > 0.0) || (xU < 0.0)
        return pow_1d(x, n, (n >= 0) ? n*xL^(n-1) : n*xU^(n-1))
    elseif iseven(n)
        return pow_even(x, n)
    end
    (xL < 0.0 < xU) && error("Undefined domain...")
    return pow_odd(x, n)
end

function ^(x::AffineEAGO{N}, n::Float64) where N
    a, b = bounds(x)
    (a <= 0.0) && error("Invalid domain...")
    if (n > 1.0) || (n < 0.0)
        return pow_1d(x, n, n*xU^(n-1))
    end
    return pow_1d(x, n, n*xL^(n-1))
end

# DONE...
function ^(x::AffineEAGO{N}, n::Number) where N
    if iszero(n) 
        return zero(x)
    elseif isone(n)
        return one(x)
    end
    return x^n
end

function inv(x::AffineEAGO{N}) where N
    a, b = bounds(x)
    (a < 0.0 < b) && error("Invalid domain...")
    if b < 0.0
        return -inv(-x)
    end
    if USE_MIN_RANGE
        p = -1/b^2
        q = -(p*(a + b)^2)/(2*a)
        Δ = -(p*(a - b)^2)/(2*a)
        return AffineEAGO(x, p, q, Δ)
    end
    p = -1/(a*b)
    q = -0.5*p*(sqrt(a) + sqrt(b))^2
    Δ = -0.5*p*(sqrt(a) - sqrt(b))^2
    return AffineEAGO(x, p, q, Δ)
end

struct MCAffPnt{Q,N,T}
    v::MC{N,T}
    box::AffineEAGO{N}
end

zero(::Type{MCAffPnt{Q,N,T}}) where {Q,N,T} = MCAffPnt{Q,N,T}(zero(MC{N,T}), zero(AffineEAGO{N}))
zero(x::MCAffPnt{Q,N,T}) where {Q,N,T} = MCAffPnt{Q,N,T}(zero(x.v), zero(x.box))

+(x::MCAffPnt{Q,N,T}, y::MCAffPnt{Q,N,T}) where {Q,N,T} = MCAffPnt{Q,N,T}(x.v + y.v, x.box + y.box)
+(x::MCAffPnt{Q,N,T}, α::T) where {Q,N,T} = MCAffPnt{Q,N,T}(x.v + α, x.box + α)
+(α::T, x::MCAffPnt{Q,N,T}) where {Q,N,T} = x + α

*(x::MCAffPnt{Q,N,T}, y::MCAffPnt{Q,N,T}) where {Q,N,T} = MCAffPnt{Q,N,T}(x.v*y.v, x.box*y.box)
*(x::MCAffPnt{Q,N,T}, α::T) where {Q,N,T} = MCAffPnt{Q,N,T}(x.v * α, x.box * α)
*(α::T, x::MCAffPnt{Q,N,T}) where {Q,N,T} = x*α

-(x::MCAffPnt{Q,N,T}) where {Q,N,T} = MCAffPnt{Q,N,T}(-x.v, -x.box)
-(x::MCAffPnt{Q,N,T}, α::T) where {Q,N,T} = x + (-α)
-(α::T, x::MCAffPnt{Q,N,T}) where {Q,N,T} = α + (-x)

/(x::MCAffPnt{Q,N,T}, α::T) where {Q,N,T} = MCAffPnt{Q,N,T}(x.v/α, x.box/α)
/(α::T, x::MCAffPnt{Q,N,T}) where {Q,N,T} = MCAffPnt{Q,N,T}(α*inv(x.v), α*inv(x.box))

^(x::MCAffPnt{Q,N,T}, n::Integer) where {Q,N,T} = MCAffPnt{Q,N,T}(x.v^n, x.box^n)
inv(x::MCAffPnt{Q,N,T}) where {Q,N,T} = MCAffPnt{Q,N,T}(inv(x.v), inv(x.box))

log(x::MCAffPnt{Q,N,T}) where {Q,N,T} = MCAffPnt{Q,N,T}(log(x.v), log(x.box))
log10(x::MCAffPnt{Q,N,T}) where {Q,N,T} = MCAffPnt{Q,N,T}(log10(x.v), log10(x.box))

exp(x::MCAffPnt{Q,N,T}) where {Q,N,T} = MCAffPnt{Q,N,T}(exp(x.v), exp(x.box))
exp10(x::MCAffPnt{Q,N,T}) where {Q,N,T} = MCAffPnt{Q,N,T}(exp10(x.v), exp10(x.box))

function extract_apriori_info(t::RelaxAA, x::AffineEAGO{N}, y::MC{N,T}) where {Q,N,T}
    p, P = t.p, t.P
    padj = (t.p - mid(t.P))/diam(t.P)
    z = sum(i -> x.γ[i]*padj[i], Val(N))
    Z = x.c + sum(i -> x.γ[i]*UNIT_INTERVAL, Val(N))
    xcv  = x.c + z - x.Δ
    xcc  = x.c + z + x.Δ
    xcvU = Z.hi - x.Δ
    xccL = Z.lo + x.Δ
    xcvg = x.γ
    xccg = x.γ
    return cv, cvU, cc, ccL, cv_grad, cc_grad
end

relax_info(s::RelaxAA, n::Int, t::T) where T = MCAffPnt{AffineEAGO{n},n,T}
f_init!(::RelaxAA, g::DAT, b::RelaxCache) = nothing