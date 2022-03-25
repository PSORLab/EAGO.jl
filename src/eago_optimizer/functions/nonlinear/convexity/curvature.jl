#=
@enum(Vexity, VEX_CONST, VEX_AFFINE, VEX_CONVEX, VEX_CONCAVE, VEX_NONDCP, VEX_UNSET)

function -(v::Vexity)
    if v == VEX_CONVEX
        return VEX_CONCAVE
    elseif v == VEX_CONCAVE
        return VEX_CONCAVE
    end
    return v
end

function +(v::Vexity, w::Vexity)
    if (v == VEX_CONST) && (w == VEX_CONST)            # constant times other vexity cases
        return v
    elseif (v == VEX_CONST) && (w == VEX_NONDCP)
        return w
    elseif (v == VEX_NONDCP) && (w == VEX_CONST)
        return v
    elseif (v == VEX_CONST)
        return w
    elseif (w == VEX_CONST)
        return v
    elseif (v == VEX_AFFINE) && (w == VEX_AFFINE)      # affine times other vexity cases
        return v
    elseif (v == VEX_AFFINE) && (w == VEX_CONVEX)
        return w
    elseif (v == VEX_CONVEX) && (w == VEX_AFFINE)
        return v
    elseif (v == VEX_AFFINE) && (w == VEX_CONCAVE)
        return w
    elseif (v == VEX_CONCAVE) && (w == VEX_AFFINE)     # remaining vexity cases
        return v
    elseif (v == VEX_CONVEX) && (w == VEX_CONVEX)
        return v
    elseif (v == VEX_CONCAVE) && (w == VEX_CONCAVE)
        return v
    elseif (v == VEX_CONCAVE) && (w == VEX_CONVEX)
        return VEX_NONDCP
    end
    return VEX_NONDCP                                  # (v == VEX_CONVEX) && (w == VEX_CONCAVE)
end

curvature(f, x::Interval) = error("curvature(f,x) for $f is not implemented.")

function curvature(d2f::Interval{T}, x::Interval{T}) where T
    if zero(T) <= d2f
        return VEX_CONCAVE
    elseif zero(T) >= d2f
        return VEX_CONVEX
    elseif isatomic(x)
        return VEX_CONST
    end
    return VEX_NONDCP
end

for F in (positive, negative, lower_bnd, upper_bnd, bnd, rad2deg, deg2rad)
    @eval curvature(::typeof($F), x::Interval) = VEX_AFFINE
end

curvature(::typeof(+), x::Interval, y::Interval, args...) = VEX_CONST
curvature(::typeof(+), x::Interval) = VEX_CONST

curvature(::typeof(-), x::Interval, y::Interval) = VEX_AFFINE
curvature(::typeof(-), x::Interval) = VEX_AFFINE

for F in (abs, abs2, exp, exp2, exp10, expm1, leaky_relu, xlogx, cosh)
    @eval curvature(::typeof($F), x::Interval) = VEX_CONVEX
end

for F in (log, log2, log10, log1p, acosh, sqrt)
    @eval curvature(::typeof($F), x::Interval) = VEX_CONCAVE
end

curvature(::typeof(sin), x::Interval) = curvature(-sin(x), x)
curvature(::typeof(cos), x::Interval) = curvature(-cos(x), x)
curvature(::typeof(tan), x::Interval) = curvature(2*tan(x)*sec(x)^2, x)

for F in (asin, acot, sinh, csch, coth, atanh, acsch, xabsx, erfinv, erfc)
    @eval function curvature(::typeof($F), x::Interval)
        if x <= zero(T)
            return VEX_CONCAVE
        elseif x >= zero(T)
            return VEX_CONVEX
        end
        return VEX_NONDCP
    end
end

for F in (acos, atan, cbrt, tanh, erf, erfcinv)
    @eval function curvature(::typeof($F), x::Interval)
        if x <= zero(T)
            return VEX_CONVEX
        elseif x >= zero(T)
            return VEX_CONCAVE
        end
        return VEX_NONDCP
    end
end

function curvature(::typeof(xexpax), x::Interval{T}, a::Interval{T}) where T
    if isatomic(a)
        inflect_pnt = -2/a.lo
        if inflect_pnt < x
            return VEX_CONVEX
        elseif inflect_pnt > x
            return VEX_CONCAVE
        end
    end
    return VEX_NONDCP
end

# TODO: HANDLE ARGMAX/ARGMIN IS CONSTANT CASE -> VEXAFFINE CASE... (need z interval bounds?)
curvature(::typeof(max), x::Interval{T}, args...) = VEX_CONVEX
curvature(::typeof(min), x::Interval{T}, args...) = VEX_CONCAVE

function curvature(::typeof(relu), x::Interval{T}) where T
    if x >= zero(T)
        return VEX_CONVEX
    end
    return VEX_CONST
end

# TODO: sech
# TODO: maxtanh, pentanh, sigmoid, bisigmoid, softsign, softplus, elu
=#