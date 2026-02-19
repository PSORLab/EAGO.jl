
const USE_MIN_RANGE = true

struct AffineEAGO{N}
    c::Float64               # Mid-point
    γ::SVector{N,Float64}    # Affine terms
    Δ::Float64               # Error term
end

function AffineEAGO(x::AffineEAGO{N}, p::Float64, q::Float64, δ::Float64) where N
    c = p * x.c + q
    γ = p .* x.γ
    Δ = p * x.Δ + δ
    AffineEAGO{N}(c, γ, δ)
end
mid(x::Interval{Float64}) = 0.5*(x.bareinterval.lo + x.bareinterval.hi)
AffineEAGO{N}(x::Float64, X::Interval{Float64}, i::Int) where N = AffineEAGO{N}(mid(X), radius(X)*seed_gradient(i, Val(N)), 0.0)

const UNIT_INTERVAL = interval(-1.0, 1.0)
interval(x::AffineEAGO{N}) where N = x.c + x.Δ*UNIT_INTERVAL + sum(y -> abs(y)*UNIT_INTERVAL, x.γ)
function bounds(x::AffineEAGO{N}) where N
    z = interval(x)
    z.bareinterval.lo, z.bareinterval.hi
end

zero(::Type{AffineEAGO{N}}) where N = AffineEAGO{N}(0.0, zero(SVector{N,Float64}), 0.0)
zero(x::AffineEAGO{N}) where N = AffineEAGO{N}(0.0, zero(SVector{N,Float64}), 0.0)

one(::Type{AffineEAGO{N}}) where N = AffineEAGO{N}(1.0, zero(SVector{N,Float64}), 0.0)
one(x::AffineEAGO{N}) where N = AffineEAGO{N}(1.0, zero(SVector{N,Float64}), 0.0)

+(x::AffineEAGO{N}, y::AffineEAGO{N}) where N = AffineEAGO{N}(x.c + y.c, x.γ .+ y.γ, x.Δ + y.Δ)
+(x::AffineEAGO{N}, α::Float64) where N = AffineEAGO{N}(α+x.c, x.γ, x.Δ)
+(α::Float64, x::AffineEAGO{N}) where N = x + α

function *(x::AffineEAGO{N}, y::AffineEAGO{N}) where N
    x0 = x.c
    y0 = y.c
    γ = SVector{N,Float64}(ntuple(i -> x0*y.γ[i] + y0*x.γ[i], Val(N)))
    Δ = abs(x0)*y.Δ + abs(y0)*x.Δ
    sx = abs(x.Δ)
    sy = abs(y.Δ)
    for i = 1:N
        sx += abs(x.γ[i])
        sy += abs(y.γ[i])
    end
    Δ += sx + sy
    return AffineEAGO{N}(x0*y0, γ, Δ)
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


function exp(x::AffineEAGO{N}) where N
    a, b = bounds(x)
    fa = exp(a)
    fb = exp(b)
    if USE_MIN_RANGE
        p = exp(a)
        q = 0.5*(fa + fb - p*(a + b))
        Δ = abs(0.5*(fb - fa - p*(b - a)))
        return AffineEAGO(x, p, q, Δ)
    end
    p = (fb - fa)/(b - a)
    ξ = log(p)
    fξ = p
    q = 0.5*(fa + fξ - p*(a + ξ))
    Δ = abs(0.5*(fξ - fa - p*(ξ - a)))
    return AffineEAGO(x, p, q, Δ)
end

function exp10(x::AffineEAGO{N}) where N
    a, b = bounds(x)
    fa = exp10(a)
    fb = exp10(b)
    if USE_MIN_RANGE
        p = log10(a)
        q = 0.5*(fa + fb - p*(a + b))
        Δ = abs(0.5*(fb - fa - p*(b - a)))
        return AffineEAGO(x, p, q, Δ)
    end
    p = (fb - fa)/(b - a)
    ξ = log10(p)
    fξ = p
    q = 0.5*(fa + fξ - p*(a + ξ))
    Δ = abs(0.5*(fξ - fa - p*(ξ - a)))
    return AffineEAGO(x, p, q, Δ)
end

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

function pow_1d(x::AffineEAGO{N}, n::Number, p) where N
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
        m = min(0.0, fa, fb)
        M = max(0.0, fa, fb)
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
    #println("ran power odd")
    # TODO: DOES THIS HANDLE a <= 0.0 <= b?
    a, b = bounds(x)
    fa = a^n
    fb = b^n
    p = (fb - fa)/(b - a)
    q = 0.0
    ξ = (p/n)^(1/(n-1))
    fξ = ξ^n
    Δ = abs(fξ - p*ξ)
    return AffineEAGO(x, p, q, Δ)
    #y = Base.power_by_squaring(x,n)
   # return y
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
    return pow_odd(x, n)
end

function ^(x::AffineEAGO{N}, n::Float64) where N
    if isinteger(n)
        return x^Int(n)
    end
    xL, xU = bounds(x)
    (xL < 0.0) && error("Invalid domain...")
    if (n > 1.0) || (n < 0.0)
        return pow_1d(x, n, n*xU^(n-1))
    end
    return pow_1d(x, n, n*xL^(n-1))
end


^(x::AffineEAGO{N}, y::AffineEAGO{N}) where N = exp(y*log(x))
function ^(x::AffineEAGO{N}, n) where N
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

struct MCAffPnt{N,T}
    v::MC{N,T}
    box::AffineEAGO{N}
end
MC(x::MCAffPnt{N,T}) where {N,T<:RelaxTag} = x.v
MC(x::MC{N,T}) where {N, T<:RelaxTag} = x

relax_info(s::RelaxAA, n::Int, t::T) where T = MCAffPnt{n,T}

zero(::Type{MCAffPnt{N,T}}) where {N,T} = MCAffPnt{N,T}(zero(MC{N,T}), zero(AffineEAGO{N}))
zero(x::MCAffPnt{N,T}) where {N,T} = MCAffPnt{N,T}(zero(x.v), zero(x.box))

one(::Type{MCAffPnt{N,T}}) where {N,T} = MCAffPnt{N,T}(one(MC{N,T}), one(AffineEAGO{N}))
one(x::MCAffPnt{N,T}) where {N,T} = MCAffPnt{N,T}(one(MC{N,T}), one(AffineEAGO{N}))

+(x::MCAffPnt{N,T}, y::MCAffPnt{N,T}) where {N,T} = MCAffPnt{N,T}(x.v + y.v, x.box + y.box)
+(x::MCAffPnt{N,T}, α::Number) where {N,T} = MCAffPnt{N,T}(x.v + α, x.box + α)
+(α::Number, x::MCAffPnt{N,T}) where {N,T} = x + α

*(x::MCAffPnt{N,T}, y::MCAffPnt{N,T}) where {N,T} = MCAffPnt{N,T}(x.v*y.v, x.box*y.box)
*(x::MCAffPnt{N,T}, α::Number) where {N,T} = MCAffPnt{N,T}(x.v * α, x.box * α)
*(α::Number, x::MCAffPnt{N,T}) where {N,T} = x*α

-(x::MCAffPnt{N,T}, y::MCAffPnt{N,T}) where {N,T} = MCAffPnt{N,T}(x.v-y.v, x.box-y.box)
-(x::MCAffPnt{N,T}) where {N,T} = MCAffPnt{N,T}(-x.v, -x.box)
-(x::MCAffPnt{N,T}, α::Number) where {N,T} = x + (-α)
-(α::Number, x::MCAffPnt{N,T}) where {N,T} = α + (-x)

/(x::MCAffPnt{N,T}, y::MCAffPnt{N,T}) where {N,T} = x*inv(y)
/(x::MCAffPnt{N,T}, α::T) where {N,T} = MCAffPnt{N,T}(x.v/α, x.box/α)
/(α::T, x::MCAffPnt{N,T}) where {N,T} = MCAffPnt{N,T}(α*inv(x.v), α*inv(x.box))

^(x::MCAffPnt{N,T}, n::Integer) where {N,T} = MCAffPnt{N,T}(x.v^n, x.box^n)
^(x::MCAffPnt{N,T}, n::Number) where {N,T} = MCAffPnt{N,T}(x.v^n, x.box^n)
^(x::MCAffPnt{N,T}, n::MCAffPnt{N,T}) where {N,T} = MCAffPnt{N,T}(x.v^n.v, x.box^n.box)

for op in (:inv, :log, :log10, :exp, :exp10)
    @eval ($op)(x::MCAffPnt{N,T}) where {N,T} = MCAffPnt{N,T}(($op)(x.v), ($op)(x.box))
end

interval(x::MCAffPnt{N,T}) where {N,T<:RelaxTag} = intersect_interval(interval(x.v), interval(x.box))

function cut(x::MCAffPnt{N,T}, z::MCAffPnt{N,T}, v::VariableValues, ϵ::Float64, s::Vector{Int}, cflag::Bool, pflag::Bool) where {N,T<:RelaxTag}
    (pflag & cflag)  && (return set_value_post(intersect(x, interval(z)), v, s, ϵ))
    (pflag & !cflag) && (return set_value_post(x, v, s, ϵ))
    (pflag & cflag)  && (return intersect(x, interval(z)))
    return x
end

function cut(x::MC{N,T}, z::MCAffPnt{N,T}, v::VariableValues, ϵ::Float64, s::Vector{Int}, cflag::Bool, pflag::Bool) where {N,T<:RelaxTag}
    (pflag & cflag)  && (return set_value_post(intersect(x, interval(z)), v, s, ϵ))
    (pflag & !cflag) && (return set_value_post(x, v, s, ϵ))
    (pflag & cflag)  && (return intersect(x, interval(z)))
    return x
end

function varset(::Type{MCAffPnt{N,T}}, i, x_cv, x_cc, l, u) where {N,T<:RelaxTag}
    v = seed_gradient(i, Val(N))
    v_Intv = interval(l, u)
    v_mc = MC{N,T}(x_cv, x_cc, v_Intv, v, v, false) 
    v_aff = AffineEAGO{N}(x_cv, v_Intv, i)
    return MCAffPnt{N,T}(v_mc, v_aff)
end

function fprop!(t::RelaxAAInfo, vt::Variable, g::DAT, b::RelaxCache{MCAffPnt{N,T},N,T}, k) where {N,T<:RelaxTag}
    i = first_index(g, k)
    x = val(b, i)
    l = lbd(b, i)
    u = ubd(b, i)
    if l == u
        b[k] = x
        b._is_num[k] = true
    else
        z = varset(MCAffPnt{N,T}, rev_sparsity(g, i, k), x, x, l, u)
        if !first_eval(t, b)
            z = intersect(z, interval(b, k))
        end
        b._info[k] = z
        b._is_num[k] = false
    end
    nothing
end

function fprop!(t::RelaxAAInfo, ex::Subexpression, g::DAT, b::RelaxCache{MCAffPnt{N,T},N,T}, k) where {N,T<:RelaxTag}
    x =  first_index(g, k)
    if subexpression_is_num(b, x)
        b[k] = subexpression_num(b, x)
        b._is_num[k] = true
    else
        b._info[k] = subexpression_info(b, x)
        b._is_num[k] = false
    end
end

for (F, f) in ((PLUS, :+), (MIN, :min), (MAX, :max), (DIV, :/), (ARH, :arh))
    @eval function fprop_2!(t::RelaxAAInfo, v::Val{$F}, g::DAT, b::RelaxCache{MCAffPnt{N,T},N,T}, k) where {N,T<:RelaxTag}
        x = child(g, 1, k)
        y = child(g, 2, k)
        if !xy_num(b, x, y)
            if xyset(b, x, y)
                z = ($f)(info(b, x), info(b, y))
            elseif xset_ynum(b, x, y)
                z = ($f)(info(b, x), num(b, y))
            else
                z = ($f)(num(b, x), info(b, y))
            end
            b._info[k] = cut(z, info(b,k), b.ic.v, b.ϵ_sg, sparsity(g, k), b.cut, false)
            b._is_num[k] = false
        else
            b[k] = ($f)(num(b, x), num(b, y))
            b._is_num[k] = true
        end
    end
end

function fprop!(t::RelaxAAInfo, v::Val{MINUS}, g::DAT, b::RelaxCache{MCAffPnt{N,T},N,T}, k) where {N,T<:RelaxTag}
    x = child(g, 1, k)
    if is_binary(g, k)
        y = child(g, 2, k)
        if !xy_num(b, x, y)
            if xyset(b, x, y)
                z = info(b, x) - info(b, y)
            elseif xset_ynum(b, x, y)
                z = info(b, x) - num(b, y)
            else
                z = num(b, x) - info(b, y)
            end
            b._info[k] = cut(z, info(b,k), b.ic.v, b.ϵ_sg, sparsity(g, k), b.cut, false)
            b._is_num[k] = false
        else
            b[k] = num(b, x) - num(b, y)
            b._is_num[k] = true
        end
    else
        if is_num(b, x)
            b[k] = -num(b, x)
            b._is_num[k] = true
        else
            b._info[k] = cut(-info(b, x), info(b, k), b.ic.v, b.ϵ_sg, sparsity(g, k), b.cut, false)
            b._is_num[k] = false
        end
    end
end

function fprop_n!(t::RelaxAAInfo, v::Val{PLUS}, g::DAT, b::RelaxCache{MCAffPnt{N,T},N,T}, k::Int) where {N,T<:RelaxTag}
    z = zero(MCAffPnt{N,T})
    znum = 0.0
    numval = true
    for i in children(g, k)
        if is_num(b, i)
            znum += num(b, i)
        else
            numval = false
            z += info(b, i)
        end
    end
    if numval
        b[k] = znum
        b._is_num[k] = true
    else
        z += znum
        b._info[k] = cut(z, info(b, k), b.ic.v, b.ϵ_sg, sparsity(g, k), b.cut, false)
        b._is_num[k] = false
    end
end

function fprop_n!(t::RelaxAAInfo, v::Val{MIN}, g::DAT, b::RelaxCache{MCAffPnt{N,T},N,T}, k::Int) where {N,T<:RelaxTag}
    z = Inf*one(MCAffPnt{N,T})
    znum = Inf
    numval = true
    for i in children(g, k)
        if is_num(b, i)
            znum = min(znum, num(b, i))
        else
            numval = false
            z = min(z, info(b, i))
        end
    end
    if numval
        b[k] = znum
        b._is_num[k] = true
    else
        z = min(z, znum)
        b._info[k] = cut(z, info(b, k), b.ic.v, b.ϵ_sg, sparsity(g, k), b.cut, false)
        b._is_num[k] = false
    end
end

function fprop_n!(t::RelaxAAInfo, v::Val{MAX}, g::DAT, b::RelaxCache{MCAffPnt{N,T},N,T}, k::Int) where {N,T<:RelaxTag}
    z = -Inf*one(MCAffPnt{N,T})
    znum = -Inf
    numval = true
    for i in children(g, k)
        if is_num(b, i)
            znum = max(znum, num(b, i))
        else
            numval = false
            z = max(z, info(b, i))
        end
    end
    if numval
        b[k] = znum
        b._is_num[k] = true
    else
        z = max(z, znum)
        b._info[k] = cut(z, info(b, k), b.ic.v, b.ϵ_sg, sparsity(g, k), b.cut, false)
        b._is_num[k] = false
    end
end

function fprop_2!(t::RelaxAAInfo, v::Val{MULT}, g::DAT, b::RelaxCache{MCAffPnt{N,T},N,T}, k::Int) where {N,T<:RelaxTag}

    x = child(g, 1, k)
    y = child(g, 2, k)

    if !xy_num(b, x, y)
        if xyset(b, x, y)
            xr = info(b, x)
            yr = info(b, y)
            xv = xr.v
            yv = yr.v
            if b.use_apriori_mul
                dp = b.dp
                dP = b.dP
                p_rel = b.p_rel
                p_diam = b.p_diam
                s = sparsity(g, 1)
                u1max, u2max, v1nmax, v2nmax = estimator_extrema(xr, yr, s, dP)
                zv = xv*yv
                wIntv = zv.Intv
                if (u1max < xv.Intv.bareinterval.hi) || (u2max < yv.Intv.bareinterval.hi)
                    u1cv, u2cv, u1cvg, u2cvg = estimator_under(xv, yv, xr, yr, s, dp, dP, p_rel, p_diam)
                    za_l = McCormick.mult_apriori_kernel(xv, yv, wIntv, u1cv, u2cv, u1max, u2max, u1cvg, u2cvg)
                    zv = intersect(zv, za_l)
                end
                if (v1nmax > -xv.Intv.bareinterval.lo) || (v2nmax > -yv.Intv.bareinterval.lo)
                    v1ccn, v2ccn, v1ccgn, v2ccgn = estimator_over(xv, yv, xr, yr, s, dp, dP, p_rel, p_diam)
                    za_u = McCormick.mult_apriori_kernel(-xv, -yv, wIntv, v1ccn, v2ccn, v1nmax, v2nmax, v1ccgn, v2ccgn)
                    zv = intersect(zv, za_u)
                end
                z = MCAffPnt{N,T}(zv, xr.box*yr.box)
            else
                z = xv*yv
            end
        elseif xset_ynum(b, x, y)
            z = info(b, x)*num(b, y)
        else
            z = num(b, x)*info(b, y)
        end
        b._info[k] = cut(z, info(b,k), b.ic.v, b.ϵ_sg, sparsity(g, k), b.cut, false)
        b._is_num[k] = true
    else
        b[k] = num(b, x)*num(b, y)
        b._is_num[k] = false
    end
end

function fprop_n!(t::RelaxAAInfo, ::Val{MULT}, g::DAT, b::RelaxCache{MCAffPnt{N,T},N,T}, k::Int) where {N,T<:RelaxTag}

    z = one(MCAffPnt{N,T})
    znum = one(Float64)
    numval = true
    if b.use_apriori_mul
        zr = one(MCAffPnt{N,T})
        dp = b.dp
        dP = b.dP
        s = sparsity(g, 1)
        for (q,i) in enumerate(children(g, k))
            if is_num(b, i)
                znum = znum*num(b, i)
            else
                numval = false
                xi = info(b, i)
                x = xi.v
                xr = info(b, i)
                u1max, u2max, v1nmax, v2nmax = estimator_extrema(zr, xr, s, dP)
                zv = z*x
                wIntv = zv.Intv
                if (u1max < z.Intv.bareinterval.hi) || (u2max < x.Intv.bareinterval.hi)
                    u1cv, u2cv, u1cvg, u2cvg = estimator_under(zr, xr, s, dp, dP)
                    za_l = McCormick.mult_apriori_kernel(z, x, wIntv, u1cv, u2cv, u1max, u2max, u1cvg, u2cvg)
                    zv = intersect(zv, za_l)
                end
                if (v1nmax > -z.Intv.bareinterval.lo) || (v2nmax > -x.Intv.bareinterval.lo)
                    v1ccn, v2ccn, v1ccgn, v2ccgn = estimator_under(zr, xr, s, dp, dP)
                    za_u = McCormick.mult_apriori_kernel(-z, -x, wIntv, v1ccn, v2ccn, v1nmax, v2nmax, v1ccgn, v2ccgn)
                    zv = intersect(zv, za_u)
                end
                zr = zr*xr
                zv = cut(zv, zv, b.ic.v, b.ϵ_sg, sparsity(g, i), b.cut, false)
                z = zv
            end
        end
    else
        for (q,i) in enumerate(children(g, k))
            if is_num(b, i)
                znum = znum*num(b, i)
            else
                numval = false
                x = info(b, i)
                z = z*x
                z = cut(z, z, b.ic.v, b.ϵ_sg, sparsity(g, i), b.cut, false)
            end
        end
    end
    if numval
        b[k] = znum
        b._is_num[k] = true
    else
        z = z*znum
        b._info[k] = z
        b._is_num[k] = false
    end
end

for F in (PLUS, MULT, MIN, MAX)
    @eval function fprop!(t::RelaxAAInfo, v::Val{$F}, g::DAT, b::RelaxCache{MCAffPnt{N,T},N,T}, k) where {N,T<:RelaxTag}
        is_binary(g, k) ? fprop_2!(RelaxAAInfo(), Val($F), g, b, k) : fprop_n!(RelaxAAInfo(), Val($F), g, b, k)
    end
end
function fprop!(t::RelaxAAInfo, v::Val{DIV}, g::DAT, b::RelaxCache{MCAffPnt{N,T},N,T}, k) where {N,T<:RelaxTag}
    fprop_2!(RelaxAAInfo(), v, g, b, k)
end

function fprop!(t::RelaxAAInfo, v::Val{POW}, g::DAT, b::RelaxCache{MCAffPnt{N,T},N,T}, k::Int) where {N,T<:RelaxTag}
    x = child(g, 1, k)
    y = child(g, 2, k)
    if is_num(b, y) && isone(num(b, y))
        b._info[k] = info(b, x)
    elseif is_num(b,y) && iszero(num(b, y))
        b._info[k] = one(MCAffPnt{N,T})
    elseif !xy_num(b, x, y)
        if xyset(b, x, y)
            z = info(b, x)^info(b, y)
        elseif xset_ynum(b, x, y)
            z = info(b, x)^num(b, y)
        else
            z = num(b, x)^info(b, y)
        end
        b._info[k] = cut(z, info(b,k), b.ic.v, b.ϵ_sg, sparsity(g, k), b.cut, false)
        b._is_num[k] = false
    else
        b[k] = num(b, x)^num(b, y)
        b._is_num[k] = true
    end
end

function fprop!(t::RelaxAAInfo, v::Val{USER}, g::DAT, b::RelaxCache{MCAffPnt{N,T},N,T}, k::Int) where {N,T<:RelaxTag}
    f = user_univariate_operator(g, first_index(g, k))
    x = child(g, 1, k)
    if is_num(b, x)
        b[k] = f(num(b, x))
        b._is_num[k] = true
    else
        z = f(info(b, x))
        b._info[k] = cut(z, info(b, k), b.ic.v, zero(Float64), sparsity(g, k), b.cut, b.post)
        b._is_num[k] = false
    end
end

function fprop!(t::RelaxAAInfo, v::Val{USERN}, g::DAT, b::RelaxCache{MCAffPnt{N,T},N,T}, k::Int) where {N,T<:RelaxTag}
    mv = user_multivariate_operator(g, first_index(g, k))
    n = arity(g, k)
    set_input = _info_input(b, n)
    num_input = _num_input(b, n)
    anysets = false
    i = 1
    for c in children(g, k)
        if is_num(b, c)
            x = num(b, c)
            if !isinf(x)
                set_input[i] = MCAffPnt{N,T}(x)
                num_input[i] = x
            end
        else
            set_input[i] = info(b, c)
            anysets = true
        end
        i += 1
    end
    if anysets
        z = MOI.eval_objective(mv, set_input)::MCAffPnt{N,T}
        b._info[k] = cut(z, info(b, k), b.ic.v, zero(Float64), sparsity(g,k), b.cut, b.post)
        b._is_num[k] = false
    else
        b[k] = MOI.eval_objective(mv, num_input)
        b._is_num[k] = true
    end
end

for ft in UNIVARIATE_ATOM_TYPES
    f = UNIVARIATE_ATOM_DICT[ft]
    (f == :user || f == :+ || f == :-) && continue
    @eval function fprop!(t::RelaxAAInfo, v::Val{$ft}, g::DAT, b::RelaxCache{MCAffPnt{N,T},N,T}, k) where {N,T<:RelaxTag}
        x = child(g, 1, k)
        if is_num(b, x)
            b._is_num[k] = true
            return b[k] = ($f)(num(b, x))
        else
            b._is_num[k] = false
            z = ($f)(info(b, x))
            b._info[k] = cut(z, info(b, k), b.ic.v, zero(Float64), sparsity(g,k), b.cut, b.post)
        end
    end
end

for (F, f) in ((LOWER_BND, :lower_bnd), (UPPER_BND, :upper_bnd))
    @eval function fprop!(t::RelaxAAInfo, v::Val{$F}, g::DAT, b::RelaxCache{MCAffPnt{N,T},N,T}, k) where {N,T<:RelaxTag}
        x = child(g, 1, k)
        y = child(g, 2, k)
        if is_num(b, y)
            z = info(b, x)
            z = ($f)(z, num(b, y))
            b._info[k] = cut(z, info(b, k), b.ic.v, zero(Float64), sparsity(g, k), b.cut, b.post)
        end
        b._is_num[k] = b._is_num[x]
    end
end

function fprop!(t::RelaxAAInfo, v::Val{BND}, g::DAT, b::RelaxCache{MCAffPnt{N,T},N,T}, k) where {N,T<:RelaxTag}
    x = child(g, 1, k)
    z = info(b, x)
    y = child(g, 2, k)
    r = child(g, 3, k)
    if is_num(b, y) && is_num(b, r)
        z = bnd(z, num(b, y),num(b, r))
    end
    b._is_num[k] = b._is_num[x]
    b._info[k] = cut(z, info(b, k), b.ic.v, zero(Float64), sparsity(g,k), b.cut, b.post)
end

function f_init!(t::RelaxAA, g::DAT, b::RelaxCache)
    tinfo = RelaxAAInfo()
    for k = node_count(g):-1:1
        c = node_class(g, k)
        (c == EXPRESSION)    && fprop!(tinfo, Expression(), g, b, k)
        (c == VARIABLE)      && fprop!(tinfo, Variable(), g, b, k)
        (c == SUBEXPRESSION) && fprop!(tinfo, Subexpression(), g, b, k)
        k_is_num = b._is_num[k]
        b[k] = b._info[k].v
        b._is_num[k] = k_is_num
    end
    fprop!(Relax(), g, b)
    nothing
end

for d in ALL_ATOM_TYPES
    @eval function fprop!(t::RelaxAA, v::Val{$d}, g::DAT, b::RelaxCache{V,N,T}, k::Int) where {V,N,T<:RelaxTag}
        fprop!(Relax(), v, g, b, k)
    end
end

function estimator_extrema(x::MCAffPnt{N,T}, y::MCAffPnt{N,T}, s, dP) where {N,T}

    xIntv = x.box.c + sum(z -> z*UNIT_INTERVAL, x.box.γ)
    xcvU = xIntv.bareinterval.hi - x.box.Δ
    xccL = xIntv.bareinterval.lo + x.box.Δ

    yIntv = y.box.c + sum(z -> z*UNIT_INTERVAL, y.box.γ)
    ycvU = yIntv.bareinterval.hi - y.box.Δ
    yccL = yIntv.bareinterval.lo + y.box.Δ

    return xcvU, ycvU, -xccL, -yccL
end

function estimator_under(xv, yv, x::MCAffPnt{N,T}, y::MCAffPnt{N,T}, s, dP, dp, p_rel, p_diam) where {N,T}
    x_cv = x.box.c - x.box.Δ
    y_cv = y.box.c - y.box.Δ
    for (i,p_i) ∈ enumerate(s)
        rp = 2.0*p_rel[p_i]
        x_cv += x.box.γ[i]*rp
        y_cv += y.box.γ[i]*rp
    end
    x_cv_grad = SVector{N,Float64}(ntuple(i -> 2.0*x.box.γ[i].*p_diam[s[i]], Val(N)))
    y_cv_grad = SVector{N,Float64}(ntuple(i -> 2.0*y.box.γ[i].*p_diam[s[i]], Val(N)))
    x_cv, y_cv, x_cv_grad, y_cv_grad
end

function estimator_over(xv, yv, x::MCAffPnt{N,T}, y::MCAffPnt{N,T}, s, dp, dP, p_rel, p_diam) where {N,T}
    x_cc = x.box.c + x.box.Δ
    y_cc = y.box.c + y.box.Δ
    for (i,p_i) ∈ enumerate(s)
        rp = 2.0*p_rel[p_i]
        x_cc += x.box.γ[i]*rp
        y_cc += x.box.γ[i]*rp
    end
    x_ccn_grad = SVector{N,Float64}(ntuple(i -> -2.0*x.box.γ[i].*p_diam[s[i]], Val(N)))
    y_ccn_grad = SVector{N,Float64}(ntuple(i -> -2.0*y.box.γ[i].*p_diam[s[i]], Val(N)))
    -x_cc, -y_cc, x_ccn_grad, y_ccn_grad
end


#=
function extract_apriori_info(t::Union{RelaxAA,RelaxAAInfo}, v::VariableValues{Float64}, x::AffineEAGO{N}) where {N,T}
    z = 0.0
    for (k,i) in enumerate(t.v)
        l = _bd(v, i)
        u = ubd(v, i)
        z += x.γ[k]*(_val(v, i) - 0.5*(l + u))/(u - l)
    end
    Z = x.c + sum(z -> z*UNIT_INTERVAL, x.γ)
    xcv  = x.c + z - x.Δ
    xcc  = x.c + z + x.Δ
    xcvU = Z.hi - x.Δ
    xccL = Z.lo + x.Δ
    xcvg = x.γ
    xccg = x.γ
    return xcv, xcvU, xcc, xccL, xcvg, xccg
end

function _cut_info(v::VariableValues{Float64}, z::MC{N,T}, x::MCAffPnt{N,T}) where {N,T} 
    xcv, xcvU, xcc, xccL, xcvg, xccg = extract_apriori_info(RelaxAA(), v, x.box)
    zaff = MC{N,T}(xcv, xcc, Interval(x.box), xcvg, xccg, false)
    return zaff ∩ z
end
function _cut_info(v::VariableValues{Float64}, z::MCAffPnt{N,T}, x::MCAffPnt{N,T}) where {N,T} 
    xcv, xcvU, xcc, xccL, xcvg, xccg = extract_apriori_info(RelaxAA(), v, x.box)
    zaff = MC{N,T}(xcv, xcc, Interval(x.box), xcvg, xccg, false)
    return zaff ∩ z.v
end

function _cut(t::RelaxAAInfo, b, k, x::MCAffPnt{N,T}, z::MCAffPnt{N,T}, v::VariableValues, ϵ::Float64, s::Vector{Int}, c::Bool, p::Bool) where {N,T<:RelaxTag}
    xt = p ? MCAffPnt{N,T}(set_value_post(x.v, v, s, ϵ), x.box) : x
    xtmc = _cut_info(v, xt, xt)
    xtaffp = MCAffPnt{N,T}(xtmc, xt.box)
    b._info[k] = xtaffp
    return
end

function _cut(t::RelaxAA, b, k, x::MC{N,T}, z::MCAffPnt{N,T}, v::VariableValues, ϵ::Float64, s::Vector{Int}, c::Bool, p::Bool) where {N,T<:RelaxTag}
    xMC = set_value_post(c ? x ∩ x.Intv ∩ Interval(z.box) : x, v, s, ϵ)
    xt = p ? xMC : x
    zt = _cut_info(v, xt, info(b, k))
    _store_set!(b, zt, k)
    return
end

function fprop!(t::RelaxAAInfo, vt::Variable, g::DAT, b::RelaxCache{V,N,T}, k::Int) where {V,N,T<:RelaxTag}
    i = first_index(g, k)
    x = val(b, i)
    l = lbd(b, i)
    u = ubd(b, i)
    z = varset(MC{N,T}, rev_sparsity(g, i, k), x, x, l, u)
    zaff = AffineEAGO{N}(x, Interval(l,u), i)
    zinfo = MCAffPnt{N,T}(z, zaff)
    b._info[k] = zinfo
    return 
end
function fprop!(t::RelaxAA, vt::Variable, g::DAT, b::RelaxCache{V,N,T}, k::Int) where {V,N,T<:RelaxTag}
    i = first_index(g, k)
    x = val(b, i)
    l = lbd(b, i)
    u = ubd(b, i)
    z = varset(MC{N,T}, rev_sparsity(g, i, k), x, x, l, u)
    z = z ∩ _interval(b, k)
    _store_set!(b, z, k)
end
=#







