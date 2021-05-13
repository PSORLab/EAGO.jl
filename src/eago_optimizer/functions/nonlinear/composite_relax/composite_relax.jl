struct Relax <: AbstractCacheAttribute end

Base.@kwdef mutable struct RelaxCache{V,S<:Real} <: AbstractCache
    v::VariableValues{S}            = VariableValues{S}()
    _set::Vector{V}                 = V[]
    _num::Vector{S}                 = S[]
    _is_num::Vector{Bool}           = Bool[]
    _subexpression_value::Vector{V} = V[]
    _cv_grad_buffer::Vector{S}       = S[]
    _cc_grad_buffer::Vector{S}       = S[]
    _set_mv_buffer::Vector{V}        = V[]
    post::Bool                      = false
    cut::Bool                       = false
    cut_interval::Bool              = false
    ϵ_sg::Bool                      = false
    first_eval::Bool                = false
    ctx::GuardCtx                   = GuardCtx()
end
function RelaxCache{V,S}(n::Int, m::Int, p::Int) where {V,S<:Real}
    RelaxCache{V,S}(_set                 = zeros(V, n),
                    _num                 = zeros(S, n),
                    _is_num              = zeros(Bool, n),
                    _subexpression_value = zeros(V, m),
                    _cv_grad_buffer      = zeros(S, p),
                    _cc_grad_buffer      = zeros(S, p),
                    _set_mv_buffer       = zeros(V, p))
end
function initialize!(c::RelaxCache{V,S}, g::DirectedTree{S}) where {V,S<:Real}
    n = _node_count(g)
    m = _dep_subexpr_count(g)
    p = length(_sparsity(g, 1))
    c._set                 = zeros(V, n)
    c._num                 = zeros(S, n)
    c._is_num              = zeros(Bool, n)
    c._subexpression_value = zeros(V, m)
    c._cv_grad_buffer      = zeros(S, p)
    c._cc_grad_buffer      = zeros(S, p)
    c._set_mv_buffer       = zeros(V, p)
    return
end

@inline _set(b::RelaxCache{V,S}, i::Int) where {V,S}      = b._set[i] #@inbounds b._set[i]
@inline _num(b::RelaxCache{V,S}, i::Int) where {V,S}      = b._num[i] #@inbounds b._num[i]
@inline _is_num(b::RelaxCache{V,S}, i::Int) where {V,S}   = b._is_num[i] #@inbounds b._is_num[i]
@inline _interval(b::RelaxCache{V,S}, i::Int) where {V,S} = Interval{S}(_set(b, i))
@inline _subexpression_value(b::RelaxCache{V,S}, i::Int) where {V,S} = b._subexpression_value[i] # @inbounds b._subexpression_value[i]

@inline function _store_set!(b::RelaxCache{V,S}, v::V, i::Int) where {V,S}
    #@inbounds b._set[i] = v
    b._set[i] = v
    return
end
@inline function _store_subexpression!(b::RelaxCache{V,S}, v::V, i::Int) where {V,S}
    #@inbounds b._subexpression_value[i] = v
    b._subexpression_value[i] = v
    return
end

@inline _first_eval(b::RelaxCache) = b.first_eval
@inline _val(b::RelaxCache{V,S}, i::Int) where {V,S} = _val(b.v, i)
@inline _lbd(b::RelaxCache{V,S}, i::Int) where {V,S} = _lbd(b.v, i)
@inline _ubd(b::RelaxCache{V,S}, i::Int) where {V,S} = _ubd(b.v, i)

@inline function _set_input(b::RelaxCache{V,S}, n::Int) where {V,S}
    return view(b.set_mv_buffer, 1:n)
end

include(joinpath(@__DIR__, "forward_propagation.jl"))
include(joinpath(@__DIR__, "reverse_propagation.jl"))
