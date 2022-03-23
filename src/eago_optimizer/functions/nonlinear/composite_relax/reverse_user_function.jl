#=
# TODO: Infer dimensionality of user function from context type
mutable struct RevRelaxMeta{T} end
Cassette.@context RelaxCtx

macro scalar_relax_rule(f, rev)
  esc(quote
    function Cassette.overdub(::RelaxCtx{RelaxMeta{T}}, ::typeof($f), x::VRev{T,F}) where {T<:Number,F}
        r = y -> x.rev(_val(x) ∩ $rev)
        return VRev{T,typeof(r)}(($f)(x), r)
    end
  end)
end

@norev_scalar(RelaxCtx, identity)
@norev_scalar(RelaxCtx, one)
@norev_scalar(RelaxCtx, zero)
@norev_scalar(RelaxCtx, transpose)

@scalar_set_rule(sqrt, y^2)
@scalar_set_rule(real, y)

@scalar_set_rule(acos, cos(y))
@scalar_set_rule(acot, cot(y))
@scalar_set_rule(acsc, csc(y))
@scalar_set_rule(asec, sec(y))
@scalar_set_rule(asin, sin(y))
@scalar_set_rule(atan, tan(y))

@scalar_set_rule(acosd, cosd(y))
@scalar_set_rule(acotd, cotd(y))
@scalar_set_rule(acscd, cscd(y))
@scalar_set_rule(asecd, secd(y))
@scalar_set_rule(asind, sind(y))
@scalar_set_rule(atand, tand(y))

@scalar_set_rule(asinh, sinh(y))
@scalar_set_rule(acosh, cosh(y))
@scalar_set_rule(atanh, tanh(y))
@scalar_set_rule(asech, sech(y))
@scalar_set_rule(acsch, csch(y))
@scalar_set_rule(acoth, coth(y))

@scalar_set_rule(sinh, asinh(y))
@scalar_set_rule(tanh, atanh(y))

# has point discontinuity (& unbounded)
@scalar_set_rule(csch, acsch(y))
@scalar_set_rule(coth, acoth(y))



function Cassette.overdub(::RelaxCtx{RelaxMeta{T}}, ::typeof(zeros), n::Int) where {T<:Number}
    return zeros(VRev(T), n)
end
function Cassette.overdub(::RelaxCtx{RelaxMeta{T}}, ::typeof(zeros), dims...) where {T<:Number}
    return zeros(VRev(T), dims...)
end

function Cassette.overdub(ctx::RelaxCtx{RelaxMeta{T}}, ::typeof(hcat), A...) where {T<:Number}
    vA = hcat(A...)
    sz = size(vA)
    vR = zeros(VRev, sz...)
    vR[:] = vA[:]
    return vR
end
function Cassette.overdub(ctx::RelaxCtx{RelaxMeta{T}}, ::typeof(vcat), A...) where {T<:Number}
    vA = vcat(A...)
    sz = size(vA)
    vR = zeros(VRev, sz...)
    vR[:] = vA[:]
    return vR
end

# Complex(), Real(), hypot, fma, muladd, rem2pi, mod, deg2rad, rad2deg, ldexp

#=

@scalar_rule cotd(x) -(π / oftype(x, 180)) * (1 + Ω ^ 2)
@scalar_rule cscd(x) -(π / oftype(x, 180)) * Ω * cotd(x)
@scalar_rule csch(x) -(coth(x)) * Ω
@scalar_rule sec(x) Ω * tan(x)
@scalar_rule secd(x) (π / oftype(x, 180)) * Ω * tand(x)
@scalar_rule sech(x) -(tanh(x)) * Ω

@scalar_rule acot(x) -(inv(1 + x ^ 2))
@scalar_rule acsc(x) -(inv(x ^ 2 * sqrt(1 - x ^ -2)))
@scalar_rule acsc(x::Real) -(inv(abs(x) * sqrt(x ^ 2 - 1)))
@scalar_rule asec(x) inv(x ^ 2 * sqrt(1 - x ^ -2))
@scalar_rule asec(x::Real) inv(abs(x) * sqrt(x ^ 2 - 1))

@scalar_rule cosd(x) -(π / oftype(x, 180)) * sind(x)
@scalar_rule cospi(x) -π * sinpi(x)
@scalar_rule sind(x) (π / oftype(x, 180)) * cosd(x)
@scalar_rule sinpi(x) π * cospi(x)
@scalar_rule tand(x) (π / oftype(x, 180)) * (1 + Ω ^ 2)

@scalar_rule sinc(x) cosc(x)

@scalar_rule round(x) zero(x)
@scalar_rule floor(x) zero(x)
@scalar_rule ceil(x) zero(x)

=#

#@scalar_set_rule(imag(x::Real) Zero()
=#