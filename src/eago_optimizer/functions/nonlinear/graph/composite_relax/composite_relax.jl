module CompositeRelax

    struct Relax <: AbstractCacheAttribute end

    Base.@kwdef mutable struct RelaxCache{V,S} <: AbstractCache
        v::VariableValues{S}
        _set::Vector{V}
        _num::Vector{S}
        _is_num::Vector{Bool}
        _subexpression_value::Vector{V}
        post::Bool                = false
        cut::Bool                 = false
        cut_interval::Bool        = false
        ϵ_sg::Bool                = false
        first_eval::Bool          = false
        ctx::GuardCtx
        cv_grad_buffer::Vector{S}
        cc_grad_buffer::Vector{S}
        set_mv_buffer::Vector{V}
    end
    RelaxCache(::T, n::Int) where T<:Real = ADCache{T}(fill(T, n), fill(T, n))

    @inline _set(b::RelaxCache{V,S}, i::Int) where {V,S}      = @inbounds b._set[i]
    @inline _num(b::RelaxCache{V,S}, i::Int) where {V,S}      = @inbounds b._num[i]
    @inline _is_num(b::RelaxCache{V,S}, i::Int) where {V,S}   = @inbounds b._is_num[i]
    @inline _interval(b::RelaxCache{V,S}, i::Int) where {V,S} = Interval{S}(_set(b, i))
    @inline _subexpression_value(b::RelaxCache{V,S}, i::Int) where {V,S} = @inbounds b._subexpression_value[i]

    @inline function _store_set!(b::RelaxCache{V,S}, v::V, i::Int) where {V,S}
        @inbounds b._set[i] = v
        return
    end
    @inline function _store_subexpression!(b::RelaxCache{V,S}, v::V, i::Int) where {V,S}
        @inbounds b._subexpression_value[i] = v
        return
    end

    @inline _first_eval(b::RelaxCache{V,S}) = first_eval
    @inline _val(b::RelaxCache{V,S}, i::Int) where {V,S} = _val(b.v, i)
    @inline _lbd(b::RelaxCache{V,S}, i::Int) where {V,S} = _lbd(b.v, i)
    @inline _ubd(b::RelaxCache{V,S}, i::Int) where {V,S} = _ubd(b.v, i)

    @inline function _set_input(b::RelaxCache{V,S}, n::Int) where {V,S}
        return view(b.set_mv_buffer, 1:n)
    end

    include(joinpath(@__DIR__, "forward_propagation.jl"))
    include(joinpath(@__DIR__, "reverse_propagation.jl"))
end
