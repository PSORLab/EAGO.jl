"""
    Relax

Used to dispatch relaxations to a standard 
"""
struct Relax <: AbstractCacheAttribute end

"""
    RelaxAA
"""
struct RelaxAA <: AbstractCacheAttribute
    v::Vector{Int} 
end
RelaxAA() = RelaxAA(Int[])

struct RelaxAAInfo <: AbstractCacheAttribute
    v::Vector{Int} 
end
RelaxAAInfo() = RelaxAAInfo(Int[])


"""
    RelaxMulEnum
"""
struct RelaxMulEnum <: AbstractCacheAttribute
    v::Vector{Int}
    use_info::Bool
end
RelaxMulEnum() = RelaxMulEnum(Int[], false)
RelaxMulEnum(x::Vector{Int64}) = RelaxMulEnum(x, false)
RelaxMulEnum(x::Bool) = RelaxMulEnum(Int[], x)

struct RelaxMulEnumInner <: AbstractCacheAttribute
    v::Vector{Int}
    use_info::Bool
end
RelaxMulEnumInner() = RelaxMulEnumInner(Int[], false)
RelaxMulEnumInner(x::Vector{Int64}) = RelaxMulEnumInner(x, false)
RelaxMulEnumInner(x::Bool) = RelaxMulEnumInner(Int[], x)


const RELAX_ATTRIBUTE = Union{Relax,RelaxAA,RelaxAAInfo,RelaxMulEnum,RelaxMulEnumInner}

mutable struct MCBoxPnt{Q,N,T}
    v::MC{N,T}
    box::Vector{MC{N,T}}
end

cv(x::MCBoxPnt{Q,N,T}) where {Q,N,T} = x.v.cv
cc(x::MCBoxPnt{Q,N,T}) where {Q,N,T} = x.v.cc
cv_grad(x::MCBoxPnt{Q,N,T}) where {Q,N,T} = x.v.cv_grad
cc_grad(x::MCBoxPnt{Q,N,T}) where {Q,N,T} = x.v.cc_grad

cv(x::MC{N,T}) where {N,T} = x.cv
cc(x::MC{N,T}) where {N,T} = x.cc
cv_grad(x::MC{N,T}) where {N,T} = x.cv_grad
cc_grad(x::MC{N,T}) where {N,T} = x.cc_grad

function zero(::Type{MCBoxPnt{Q,N,T}}) where {Q,N,T}
    MCBoxPnt{Q,N,T}(zero(MC{N,T}), zeros(MC{N,T}, Q))
end
function setindex!(d::MCBoxPnt{Q,N,T}, x::MC{N,T}, i::Int) where {Q,N,T}
    d.box[i] = x
end

Base.@kwdef mutable struct RelaxCache{V,N,T<:RelaxTag} <: AbstractCache
    v::VariableValues{Float64}             = VariableValues{Float64}()
    dp::Vector{Float64}                    = Float64[]
    dP::Vector{Interval{Float64}}          = Interval{Float64}[]
    _set::Vector{MC{N,T}}                  = MC{N,T}[]
    _num::Vector{Float64}                  = Float64[]
    _is_num::Vector{Bool}                  = Bool[]
    _info::Vector{V}                       = V[]
    _subexpression_set::Dict{Int,MC{N,T}}  = Dict{Int,MC{N,T}}()
    _subexpression_num::Dict{Int,Float64}  = Dict{Int,Float64}()
    _subexpression_is_num::Dict{Int,Bool}  = Dict{Int,Bool}()
    _subexpression_info::Dict{Int,MC{N,T}} = Dict{Int,MC{N,T}}()
    _cv_grad_buffer::Vector{Float64}       = Float64[]
    _cc_grad_buffer::Vector{Float64}       = Float64[]
    _set_mv_buffer::Vector{MC{N,T}}        = MC{N,T}[]
    _num_mv_buffer::Vector{Float64}        = Float64[]
    _info_mv_buffer::Vector{V}             = V[]
    _mult_temp::V                          = zero(V)
    ϵ_sg::Float64                          = 1E-11
    post::Bool                             = false
    cut::Bool                              = false
    cut_interval::Bool                     = false
    first_eval::Bool                       = true
    use_apriori_mul::Bool                  = false
end
function RelaxCache{N,T}(::Relax, n::Int, m::Int, p::Int) where {N,T<:RelaxTag}
    RelaxCache{N,T}(_set                 = zeros(MC{N,T}, n),
                    _num                 = zeros(Float64, n),
                    _is_num              = zeros(Bool, n),
                    _info                = zeros(MC{N,T}, n),
                    _cv_grad_buffer      = zeros(Float64, p),
                    _cc_grad_buffer      = zeros(Float64, p),
                    _set_mv_buffer       = zeros(MC{N,T}, p),
                    _num_mv_buffer       = zeros(p),
                    _info_mv_buffer      = zeros(MC{N,T}, p),
                    )
end
function initialize!(c::RelaxCache{V,N,T}, g::DirectedTree) where {V,N,T<:RelaxTag}

    n = _node_count(g)
    m = _dep_subexpr_count(g)
    p = length(_sparsity(g, 1))

    cnst = ConstantCache(n)
    initialize!(cnst, g)
    c._num                 = copy(cnst._num)
    c._is_num              = copy(cnst._is_num)
    c._set                 = zeros(MC{N,T}, n)
    c._cv_grad_buffer      = zeros(p)
    c._cc_grad_buffer      = zeros(p)
    c._set_mv_buffer       = zeros(MC{N,T}, p)
    c._num_mv_buffer       = zeros(p)
    c._mult_temp           = zero(V)

    for i = 1:n
        push!(c._info, zero(V))
    end
    for i = 1:p
        push!(c._info_mv_buffer, zero(V))
    end
    return
end

###
###
### Access functions for RelaxCache.
###
###
_set(b::RelaxCache{V,N,T}, i::Int) where {V,N,T<:RelaxTag}  = b._set[i]
_num(b::RelaxCache{V,N,T}, i::Int) where {V,N,T<:RelaxTag}  = b._num[i]
_set_or_num(b::RelaxCache{V,N,T}, i::Int) where {V,N,T<:RelaxTag} = !_is_num(b, i) ? _set(b,i) : _num(b,i) 
_info(b::RelaxCache{V,N,T}, i::Int) where {V,N,T<:RelaxTag} = b._info[i]
_is_num(b::RelaxCache{V,N,T}) where {V,N,T<:RelaxTag}       = b._is_num
_is_num(b::RelaxCache{V,N,T}, i::Int) where {V,N,T<:RelaxTag} = b._is_num[i]
_is_unlocked(b::RelaxCache{V,N,T}, i::Int) where {V,N,T<:RelaxTag} = !_is_num(b,i)
_interval(b::RelaxCache{V,N,T}, i::Int) where {V,N,T<:RelaxTag} = Interval{Float64}(_set(b, i))
_subexpression_set(b::RelaxCache{V,N,T}, i::Int) where {V,N,T<:RelaxTag} = b._subexpression_set[i]
_subexpression_num(b::RelaxCache{V,N,T}, i::Int) where {V,N,T<:RelaxTag} = b._subexpression_num[i]
_subexpression_is_num(b::RelaxCache{V,N,T}, i::Int) where {V,N,T<:RelaxTag} = b._subexpression_is_num[i]
function _store_num!(b::RelaxCache{V,N,T}, v::Float64, i::Int) where {V,N,T<:RelaxTag}
    b._is_num[i] = true
    b._num[i] = v
    return
end
function _store_set!(b::RelaxCache{V,N,T}, v::MC{N,T}, i::Int) where {V,N,T<:RelaxTag}
    b._is_num[i] = false
    b._set[i] = v
    return
end
_store_info!(b::RelaxCache{V,N,T}, v::V, i::Int) where {V,N,T<:RelaxTag} = (b._info[i] = v; return)

function _store_subexpression_set!(b::RelaxCache{V,N,T}, v::MC{N,T}, i::Int) where {V,N,T<:RelaxTag}
    b._subexpression_is_num[i] = false
    b._subexpression_set[i] = v
    return
end
function _store_subexpression_num!(b::RelaxCache{V,N,T}, v::MC{N,T}, i::Int) where {V,N,T<:RelaxTag}
    b._subexpression_is_num[i] = true
    b._subexpression_num[i] = v
    return
end

_first_eval(b::RelaxCache) = b.first_eval
_val(b::RelaxCache{V,N,T}, i::Int) where {V,N,T<:RelaxTag} = _val(b.v, i)
_lbd(b::RelaxCache{V,N,T}, i::Int) where {V,N,T<:RelaxTag} = _lbd(b.v, i)
_ubd(b::RelaxCache{V,N,T}, i::Int) where {V,N,T<:RelaxTag} = _ubd(b.v, i)
_set_input(b::RelaxCache{V,N,T}, n::Int) where {V,N,T<:RelaxTag} = view(b._set_mv_buffer, 1:n)
_num_input(b::RelaxCache{V,N,T}, n::Int) where {V,N,T<:RelaxTag} = view(b._num_mv_buffer, 1:n)

include(joinpath(@__DIR__, "utilities.jl"))
include(joinpath(@__DIR__, "forward_propagation.jl"))
include(joinpath(@__DIR__, "reverse_propagation.jl"))

function fprop!(t::RELAX_ATTRIBUTE, g::DAT, b::RelaxCache{V,N,T}) where {V,N,T<:RelaxTag}
    for k = _node_count(g):-1:1
        if _is_unlocked(b, k)
            c = _node_class(g, k)
            (c == EXPRESSION)    && (fprop!(t, Expression(), g, b, k);    continue)
            (c == VARIABLE)      && (fprop!(t, Variable(), g, b, k);      continue)
            (c == SUBEXPRESSION) && (fprop!(t, Subexpression(), g, b, k); continue)
        end
    end
    return
end

function rprop!(t::RELAX_ATTRIBUTE, g::DAT, b::RelaxCache{V,N,T}) where {V,N,T<:RelaxTag}
    flag = r_init!(t, g, b)
    for k = 1:_node_count(g)
        if _is_unlocked(b, k)
            nt = _node_class(g, k)
            if nt === EXPRESSION
                flag = rprop!(t, Expression(), g, b, k)
            elseif nt === VARIABLE
                flag = rprop!(t, Variable(), g, b, k)
            elseif nt === SUBEXPRESSION
                flag = rprop!(t, Subexpression(), g, b, k)
            end
        end
    end
    return flag
end