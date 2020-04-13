Cassette.@context GuardCtx

struct GuardTracker
    domain_tol::Float64
end

function Cassette.overdub(ctx::GuardCtx, ::typeof(/), x::MC{N,T}, y::MC{N,T}) where {N,T<:RelaxTag}
    if (y.Intv.lo <= ctx.metadata.domain_tol) && (y.Intv.hi <= ctx.metadata.domain_tol)
        z = union(x/MC{N,T}(Interval{Float64}(y.Intv.lo, -ctx.metadata.domain_tol)),
                  x/MC{N,T}(Interval{Float64}(ctx.metadata.domain_tol, y.Intv.hi)))
    else
        z = x/y
    end
    z
end

function Cassette.overdub(ctx::GuardCtx, ::typeof(^), x::MC{N,T}, y::MC{N,T}) where {N,T<:RelaxTag}
    x^y
end
function Cassette.overdub(ctx::GuardCtx, ::typeof(^), x::MC{N,T}, y::Float64) where {N,T<:RelaxTag}
    x^y
end
function Cassette.overdub(ctx::GuardCtx, ::typeof(^), x::Float64, y::MC{N,T}) where {N,T<:RelaxTag}
    x^y
end

for f in (log, log2, log10, sqrt)
    @eval function Cassette.overdub(ctx::GuardCtx, ::typeof($f), x::MC{N,T}) where {N,T}
        if (x.Intv.lo <= ctx.metadata.domain_tol)
            z = ($f)(MC{N,T}(Interval{Float64}(ctx.metadata.domain_tol, x.Intv.hi)))
        else
            z = ($f)(x)
        end
        z
    end
end

for f in (log1p, acosh)
    @eval function Cassette.overdub(ctx::GuardCtx, ::typeof($f), x::MC{N,T}) where {N,T<:RelaxTag}
        if (x.Intv.lo <= -1.0 + ctx.metadata.domain_tol)
            z = ($f)(MC{N,T}(Interval{Float64}(-1.0 + ctx.metadata.domain_tol, x.Intv.hi)))
        else
            z = ($f)(x)
        end
        z
    end
end

for f in (acos, asin, atanh)
    @eval function Cassette.overdub(ctx::GuardCtx, ::typeof($f), x::MC{N,T}) where {N,T<:RelaxTag}
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
