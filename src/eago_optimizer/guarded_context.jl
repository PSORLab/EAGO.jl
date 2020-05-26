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
end

function Cassette.overdub(ctx::GuardCtx, ::typeof(/), x::MC{N,T}, y::MC{N,T}) where {N, T<:RelaxTag}
    if (y.Intv.lo <= -ctx.metadata.domain_tol) && (y.Intv.hi >= ctx.metadata.domain_tol)
        z = MC{N,T}(union(x.Intv/Interval{Float64}(y.Intv.lo, -ctx.metadata.domain_tol),
                          x.Intv/Interval{Float64}(ctx.metadata.domain_tol, y.Intv.hi)))
    else
        z = x/y
    end
    z
end

function Cassette.overdub(ctx::GuardCtx, ::typeof(/), x::Float64, y::MC{N,T}) where {N, T<:RelaxTag}
    if (y.Intv.lo <= -ctx.metadata.domain_tol) && (y.Intv.hi >= ctx.metadata.domain_tol)
        z = MC{N,T}(union(x.Intv/Interval{Float64}(y.Intv.lo, -ctx.metadata.domain_tol),
                          x.Intv/Interval{Float64}(ctx.metadata.domain_tol, y.Intv.hi)))
    else
        z = x/y
    end
    z
end

function Cassette.overdub(ctx::GuardCtx, ::typeof(^), x::MC{N,T}, y::MC{N,T}) where {N, T<:RelaxTag}
    if (y < 0.0) && ((x.Intv.lo <= -ctx.metadata.domain_tol) && (x.Intv.hi >= ctx.metadata.domain_tol))
        z = MC{N,T}(union(Interval{Float64}(x.Intv.lo, -ctx.metadata.domain_tol),
                          Interval{Float64}(ctx.metadata.domain_tol, x.Intv.hi))^y)
    else
        z = x/y
    end
    z
end
function Cassette.overdub(ctx::GuardCtx, ::typeof(^), x::MC{N,T}, y::Float64) where {N, T<:RelaxTag}
    if (y < 0.0) && ((x.Intv.lo <= -ctx.metadata.domain_tol) && (x.Intv.hi >= ctx.metadata.domain_tol))
        z = MC{N,T}(union(Interval{Float64}(x.Intv.lo, -ctx.metadata.domain_tol),
                          Interval{Float64}(ctx.metadata.domain_tol, x.Intv.hi))^y)
    else
        z = x/y
    end
    z
end

for f in (log, log2, log10, sqrt)
    @eval function Cassette.overdub(ctx::GuardCtx, ::typeof($f), x::MC{N,T}) where {N, T<:RelaxTag}
        if (x.Intv.lo <= ctx.metadata.domain_tol)
            z = ($f)(MC{N,T}(Interval{Float64}(ctx.metadata.domain_tol, x.Intv.hi)))
        else
            z = ($f)(x)
        end
        z
    end
end

for f in (log1p, acosh)
    @eval function Cassette.overdub(ctx::GuardCtx, ::typeof($f), x::MC{N,T}) where {N, T<:RelaxTag}
        if (x.Intv.lo <= -1.0 + ctx.metadata.domain_tol)
            z = ($f)(MC{N,T}(Interval{Float64}(-1.0 + ctx.metadata.domain_tol, x.Intv.hi)))
        else
            z = ($f)(x)
        end
        z
    end
end

for f in (acos, asin, atanh)
    @eval function Cassette.overdub(ctx::GuardCtx, ::typeof($f), x::MC{N,T}) where {N, T<:RelaxTag}
        if (x.Intv.lo <= ctx.metadata.domain_tol - 1.0) ||
           (x.Intv.hi >= 1.0 - ctx.metadata.domain_tol)
            z = ($f)(MC{N,T}(Interval{Float64}(max(y.Intv.lo, ctx.metadata.domain_tol - 1.0),
                                               min(1.0 - ctx.metadata.domain_tol, y.Intv.hi))))
        else
            z = ($f)(x)
        end
        z
    end
end
