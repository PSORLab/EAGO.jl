
# Defines the Vrev object which holds a value and a reverse function
"""
    Vrev{T<:Number,F}

Structure used to wrap the value and reverse function
"""
struct Vrev{T<:Number,F}
  val::T
  rev::F
end
Vrev(x::T) where T = Vrev{T,typeof(identity)}(x, identity)
@inline _val(d::Vrev{T,F}) where {T,F} = d.val

zero(::Vrev{T,F}) where {T<:Number,F} = Vrev(zero(T))

macro norev_scalar(c, f)
  esc(quote
    function Cassette.overdub(::$c, ::typeof($f), x::Vrev{T,F}) where {T<:Number,F}
        r = y -> x.rev(_val(x))
        return Vrev{T,typeof(r)}(($f)(x), r)
    end
  end)
end

macro no_prop(f)
  esc(quote
    function Cassette.overdub(::$c, ::typeof($f), x::Vrev{T,F}) where {T<:Number,F}
        r = y -> x.rev(_val(x))
        return Vrev{T,typeof(r)}(($f)(x), r)
    end
  end)
end
