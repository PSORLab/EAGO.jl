# define primitives for associative terms
Cassette.overdub(ctx::TraceCtx, ::typeof(*), x, y) = afoldl(x, y)
Cassette.overdub(ctx::TraceCtx, ::typeof(afoldl), x, y, z) = afoldl(x, y, z)
Cassette.overdub(ctx::TraceCtx, ::typeof(afoldl), a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p...) = afoldl(a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p...)

# create arrays
#function Cassette.overdub(ctx::TraceCtx, ::typeof(fill!), A::Array, val)
#    fill!(EAGO.Tracer.SetTrace[SetTrace(71963392), SetTrace(72031760)], EAGO.Tracer.SetTrace(0))
#end

# primitive for array access
Cassette.overdub(ctx::TraceCtx, ::typeof(getindex), A::Array, i::Int) = getindex(A,i)
Cassette.overdub(ctx::TraceCtx, ::typeof(getindex), A::SetTraceSto, i::Int) = getindex(A,i)

# primitive for setting index of storage array
#function Cassette.overdub(ctx::TraceCtx, ::typeof(setindex!), A::Array, val::Any, i::Int)
#
#end

# handles type assertation on right-hand side for SetTrace
function Cassette.overdub(ctx::TraceCtx, ::typeof(typeassert), x::Real, type::Type)
    if ~isa(x,SetTrace)
        typeassert(x,type)
    end
    return x
end

function Cassette.prehook(ctx::TraceCtx, f, args...)
    #println(f, args)
end
