# Copyright (c) 2018: Matthew Wilhelm & Matthew Stuber.
# This work is licensed under the Creative Commons Attribution-NonCommercial-
# ShareAlike 4.0 International License. To view a copy of this license, visit
# http://creativecommons.org/licenses/by-nc-sa/4.0/ or send a letter to Creative
# Commons, PO Box 1866, Mountain View, CA 94042, USA.
#############################################################################
# EAGO
# A development environment for robust and global optimization
# See https://github.com/PSORLab/EAGO.jl
#############################################################################
# src/eago_optimizer/guarded_context.jl
# Provides utilities for dealing with nonlinear expressions that contain domain
# violations. The assumption is that domain violations only occur due to
# expansiveness of the bounds of the nonlinear terms not the underlying model.
#############################################################################

Cassette.@context GuardCtx

struct GuardTracker
    domain_tol::Float64
    guard_on::Bool
end

const IntFltIntv = Union{Int16, Int32, Int64, Float16, Float32, Float64, Interval{Float64}}

for f in (+, *, -, max, min)
    @eval Cassette.overdub(ctx::GuardCtx, ::typeof($f), x::MC{N,T}, y::MC{N,T}) where {N, T<:RelaxTag} = f(x, y)
    @eval Cassette.overdub(ctx::GuardCtx, ::typeof($f), x::S, y::S) where {S <: IntFltIntv} =  f(x,y)
    @eval Cassette.overdub(ctx::GuardCtx, ::typeof($f), x::S, y::MC{N,T}) where {N, T<:RelaxTag, S<:IntFltIntv} = f(x,y)
    @eval Cassette.overdub(ctx::GuardCtx, ::typeof($f), x::MC{N,T}, y::S) where {N, T<:RelaxTag, S<:IntFltIntv} = f(x,y)
end

function Cassette.overdub(ctx::GuardCtx, ::typeof(/), x::MC{N,T}, y::MC{N,T}) where {N, T<:RelaxTag}

    m = ctx.metadata
    if m.guard_on && (y.Intv.lo <= -m.domain_tol) && (y.Intv.hi >= m.domain_tol)
        z = MC{N,T}(union(x.Intv/Interval{Float64}(y.Intv.lo, -m.domain_tol),
                          x.Intv/Interval{Float64}(m.domain_tol, y.Intv.hi)))
    else
        z = x/y
    end

    return z
end

function Cassette.overdub(ctx::GuardCtx, ::typeof(/), x::Float64, y::MC{N,T}) where {N, T<:RelaxTag}

    m = ctx.metadata
    if m.guard_on && (y.Intv.lo <= -m.domain_tol) && (y.Intv.hi >= m.domain_tol)
        z = MC{N,T}(union(x.Intv/Interval{Float64}(y.Intv.lo, -m.domain_tol),
                          x.Intv/Interval{Float64}(m.domain_tol, y.Intv.hi)))
    else
        z = x/y
    end

    return z
end

Cassette.overdub(ctx::GuardCtx, ::typeof(/), x::S, y::S) where {S <: IntFltIntv} =  f(x,y)
Cassette.overdub(ctx::GuardCtx, ::typeof(/), x::S, y::MC{N,T}) where {N, T<:RelaxTag, S <: IntFltIntv} =  f(x,y)
Cassette.overdub(ctx::GuardCtx, ::typeof(/), x::MC{N,T}, y::S) where {N, T<:RelaxTag, S <: IntFltIntv} =  f(x,y)

function Cassette.overdub(ctx::GuardCtx, ::typeof(^), x::MC{N,T}, y::MC{N,T}) where {N, T<:RelaxTag}

    m = ctx.metadata
    if m.guard_on && (y < 0.0) && ((x.Intv.lo <= -m.domain_tol) && (x.Intv.hi >= m.domain_tol))
        z = MC{N,T}(union(Interval{Float64}(x.Intv.lo, -m.domain_tol),
                          Interval{Float64}(m.domain_tol, x.Intv.hi))^y)
    else
        z = x/y
    end

    return z
end

function Cassette.overdub(ctx::GuardCtx, ::typeof(^), x::MC{N,T}, y::Float64) where {N, T<:RelaxTag}

    m = ctx.metadata
    if m.guard_on && (y < 0.0) && ((x.Intv.lo <= -m.domain_tol) && (x.Intv.hi >= m.domain_tol))
        z = MC{N,T}(union(Interval{Float64}(x.Intv.lo, -m.domain_tol),
                          Interval{Float64}(m.domain_tol, x.Intv.hi))^y)
    else
        z = x/y
    end

    return z
end

Cassette.overdub(ctx::GuardCtx, ::typeof(^), x::S, y::S) where {S<:IntFltIntv} = f(x,y)

for f in (log, log2, log10, sqrt)
    @eval function Cassette.overdub(ctx::GuardCtx, ::typeof($f), x::MC{N,T}) where {N, T<:RelaxTag}

        m = ctx.metadata
        if m.guard_on && (x.Intv.lo <= m.domain_tol)
            z = ($f)(MC{N,T}(Interval{Float64}(m.domain_tol, x.Intv.hi)))
        else
            z = ($f)(x)
        end

        return z
    end
end

for f in (log1p, acosh)
    @eval function Cassette.overdub(ctx::GuardCtx, ::typeof($f), x::MC{N,T}) where {N, T<:RelaxTag}

        m = ctx.metadata
        if m.guard_on && (x.Intv.lo <= -1.0 + m.domain_tol)
            z = ($f)(MC{N,T}(Interval{Float64}(-1.0 + m.domain_tol, x.Intv.hi)))
        else
            z = ($f)(x)
        end

        return z
    end
end

for f in (acos, asin, atanh)
    @eval function Cassette.overdub(ctx::GuardCtx, ::typeof($f), x::MC{N,T}) where {N, T<:RelaxTag}

        m = ctx.metadata
        if m.guard_on && ((x.Intv.lo <= m.domain_tol - 1.0) ||
           (x.Intv.hi >= 1.0 - m.domain_tol))
            z = ($f)(MC{N,T}(Interval{Float64}(max(y.Intv.lo, m.domain_tol - 1.0),
                                               min(1.0 - m.domain_tol, y.Intv.hi))))
        else
            z = ($f)(x)
        end

        return z
    end
end

for f in (log, log2, log10, sqrt, log1p, acosh, acos, asin, atanh, acosd, asind)
    @eval Cassette.overdub(ctx::GuardCtx, ::typeof($f), x::S) where {S<:IntFltIntv} = f(x)
end

for f in (abs, exp, exp2, exp10, sin, tan, cos, tan, sec, csc,
          sech, csch, coth, acsch, acoth, asech, step, sign,
          asinh, tanh, atan, cosh, sind, cosd, tand, secd, cscd, cotd,
          atand, asecd, acscd, acotd, isone, isnan, empty
          convert, in, isempty, one, zero, real, eps, rad2deg, deg2rad)
    @eval Cassette.overdub(ctx::GuardCtx, ::typeof($f), x::S) where {S<:IntFltIntv} = f(x)
    @eval Cassette.overdub(ctx::GuardCtx, ::typeof($f), x::MC{N,T}) where {N,T<:RelaxTag} = f(x)
end
