@enum(Monotonicity, MONO_CONST, MONO_NONINCR, MONO_NONDECR, MONO_NONE, MONO_UNSET)

function *(m::Monotonicity, v::Vexity)
    if m == MONO_NONDECR
        return v
    elseif m == MONO_NONINCR
        return -v
    elseif m == MONO_NONE
        if v == VEX_CONVEX
            return VEX_NONDCP
        elseif v == VEX_CONCAVE
            return VEX_NONDCP
        else
            return v
        end
    end
    return VEX_NONDCP # TODO CHECK FALLBACK RETURN
end

monotonicity(f, x::Interval) = error("monotonicity(f,x) for $f is not implemented.")

function monotonicity(df::Interval{T}, x::Interval{T}) where T
    if zero(T) <= df
        return (MONO_NONDECR,)
    elseif zero(T) >= df
        return (MONO_NONINCR,)
    end
    return (MONO_NONE,)
end

for F in (positive, negative, lower_bnd, upper_bnd, bnd)
    @eval monotonicity(::typeof($F), x::Interval) = MONO_AFFINE
end

function monotonicity(::typeof(+), x::Interval, y::Interval, args...)
    ntuple(i -> MONO_NONDECR, length(args) + 2)
end
monotonicity(::typeof(+), x::Interval) = (MONO_NONDECR)

monotonicity(::typeof(-), x::Interval, y::Interval) = (MONO_NONDECR, MONO_NONINCR)
monotonicity(::typeof(-), x::Interval) = (MONO_NONINCR,)

function monotonicity(::typeof(/), x::Interval{T}, y::Interval{T}) where T
    if v > zero(T)
        return (MONO_AFFINE, MONO_NONINCR)
    elseif v < zero(T)
        return (MONO_AFFINE, MONO_NONDECR)
    end
    return (MONO_AFFINE, MONO_NONE)
end

function monotonicity(::typeof(*), args...)
    return ntuple(i -> prod_exclude(_sign, args, i)*MONO_NONDECR, length(args))
end

for F in (abs, abs2, cosh, sech)
    @eval function monotonicity(::typeof($F), x::Interval{T}) where T
        if zero(T) <= x
            return (MONO_NONDECR,)
        elseif zero(T) >= x
            return (MONO_NONINCR,)
        end
        return (MONO_NONE,)
    end
end

for F in (exp, exp2, exp10, expm1, sinh, tanh,
          log, log2, log10, log1p,
          asin, atan,
          erf, erfinv,
          sqrt, cbrt, xabsx,
          # rectifer & sigmoid activation functions
          relu, leaky_relu, maxtanh, maxsig, pentanh, sigmoid,
          bisigmoid, softsign, softplus, elu)
    @eval monotonicity(::typeof($F), x::Interval) = (MONO_NONDECR,)
end

for F in (acos, acot, erfc, erfcinv)
    @eval monotonicity(::typeof($F), x::Interval) = (MONO_NONINCR,)
end

for F in (coth, csch)
    @eval function monotonicity(::typeof($F), x::Interval{T}) where T
        if zero(T) <= x
            return (MONO_NONINCR,)
        elseif zero(T) >= x
            return (MONO_NONDECR,)
        end
        return (MONO_NONE,)
    end
end

monotonicity(::typeof(sin), x::Interval) = monotonicity(cos(x), x)
monotonicity(::typeof(cos), x::Interval) = monotonicity(sin(x), x)
monotonicity(::typeof(tan), x::Interval) = monotonicity(sec(x)^2, x)

function monotonicity(::typeof(asec), x::Interval{T}) where T
    if one(T) <= x
        return (MONO_NONDECR,)
    elseif -one(T) >= x
        return (MONO_NONINCR,)
    end
    return (MONO_NONE,)
end

function monotonicity(::typeof(acsc), x::Interval{T}) where T
    if one(T) <= x
        return (MONO_NONINCR,)
    elseif -one(T) >= x
        return (MONO_NONDECR,)
    end
    return (MONO_NONE,)
end

function monotonicity(::typeof(xlogx), x::Interval{T}) where T
    if one(T)/MathConstants.e <= x
        return (MONO_NONDECR,)
    elseif one(T)/MathConstants.e >= x
        return (MONO_NONINCR,)
    end
    return (MONO_NONE,)
end

function monotonicity(::typeof(xexpax), x::Interval{T}, a::Interval{T}) where T
    # TODO:
    return (MONO_NONE,)
end

# arity 1 implies EXPRESSION OR SELECT
function monotonicity(::Val{1}, n::NodeType, atype::AtomType, bnds::Interval)
    if n == SELECT
        error("Not implemented yet...")
    end
    if n == EXPRESSION
        # binary_switch on univariate atype
    end
end

# arity >= 2 implies EXPRESSION
function monotonicity(n::NodeType, atype::AtomType, bnds::Interval, i::Int)
    # binary_switch... atype
end

function populate_monotonicity!(n, buffer, indx)
    mono = _mono(buffer, indx)
    if _arity(n) == 1
        @inbounds mono[1] = monotonicity(Val(1),
                                         _node_type(n),
                                         _atom_type(n),
                                         _bnd(buffer, _child_id(buffer, indx)))
    else
        unsafe_map!(monotonicity, mono, n, _arity(n))
    end
    return nothing
end

#=
Gets monotonicity of buffer at i. Recalling from a cached value if
discovered or locked. Otherwise, populating the mono[i] field.
=#
function monotonicity(buffer::ConvexityBuffer, indx)
    if _discovered(buffer, indx)
        return _mono(buffer, indx)
    end
    populate_monotonicity!(_node(dag, indx), buffer, indx)
    return _mono(buffer, i)
end
