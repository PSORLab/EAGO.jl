struct Relax <: AbstractCacheAttribute end
struct RelaxAA <: AbstractCacheAttribute end
struct RelaxMulEnum <: AbstractCacheAttribute end

const RELAX_ATTRIBUTE = Union{Relax,RelaxAA,RelaxMulEnum}

Base.@kwdef mutable struct RelaxCache{V,N,T<:RelaxTag} <: AbstractCache
    v::VariableValues{Float64}             = VariableValues{Float64}()
    _set::Vector{MC{N,T}}                  = MC{N,T}[]
    _num::Vector{Float64}                  = Float64[]
    _is_num::Vector{Bool}                  = Bool[]
    _info::Vector{V}                       = V[]
    _subexpression_set::Vector{MC{N,T}}    = MC{N,T}[]
    _subexpression_info::Vector{V}         = V[]
    _cv_grad_buffer::Vector{Float64}       = Float64[]
    _cc_grad_buffer::Vector{Float64}       = Float64[]
    _set_mv_buffer::Vector{MC{N,T}}        = MC{N,T}[]
    _info_mv_buffer::Vector{V}             = V[]
    ϵ_sg::Float64                      = 1E-11
    post::Bool                         = false
    cut::Bool                          = false
    cut_interval::Bool                 = false
    first_eval::Bool                   = true
    use_apriori_mul::Bool              = false
end
function RelaxCache{N,T}(::Relax, n::Int, m::Int, p::Int) where {N,T<:RelaxTag}
    RelaxCache{N,T}(_set                 = zeros(MC{N,T}, n),
                    _num                 = zeros(Float64, n),
                    _is_num              = zeros(Bool, n),
                    _info                = zeros(MC{N,T}, n),
                    _subexpression_set   = zeros(MC{N,T}, m),
                    _subexpression_info  = zeros(MC{N,T}, m),
                    _cv_grad_buffer      = zeros(Float64, p),
                    _cc_grad_buffer      = zeros(Float64, p),
                    _set_mv_buffer       = zeros(MC{N,T}, p),
                    _info_mv_buffer      = zeros(MC{N,T}, p)
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
    c._info                = zeros(V, n)
    c._subexpression_set   = zeros(MC{N,T}, m)
    c._subexpression_info  = zeros(V, m)
    c._cv_grad_buffer      = zeros(Float64, p)
    c._cc_grad_buffer      = zeros(Float64, p)
    c._set_mv_buffer       = zeros(MC{N,T}, p)
    c._info_mv_buffer      = zeros(V, p)
    return
end

_set(b::RelaxCache{V,N,T}, i::Int) where {V,N,T<:RelaxTag}  = b._set[i]
_info(b::RelaxCache{V,N,T}, i::Int) where {V,N,T<:RelaxTag} = b._info[i]
_num(b::RelaxCache{V,N,T}, i::Int) where {V,N,T<:RelaxTag} = b._num[i] #@inbounds b._num[i]

@propagate_inbounds _is_num(b::RelaxCache{V,N,T}) where {V,N,T<:RelaxTag}         = b._is_num #@inbounds b._is_num[i]
@propagate_inbounds _is_num(b::RelaxCache{V,N,T}, i::Int) where {V,N,T<:RelaxTag} = b._is_num[i] #@inbounds b._is_num[i]

_set_or_num(b::RelaxCache{V,N,T}, i::Int) where {V,N,T<:RelaxTag} = !_is_num(b, i) ? _set(b,i) : _num(b,i) 

@propagate_inbounds _interval(b::RelaxCache{V,N,T}, i::Int) where {V,N,T<:RelaxTag} = Interval{Float64}(_set(b, i))
@propagate_inbounds _subexpression_value(b::RelaxCache{V,N,T}, i::Int) where {V,N,T<:RelaxTag} = b._subexpression_value[i] # @inbounds b._subexpression_value[i]

@propagate_inbounds _is_unlocked(b::RelaxCache, i::Int) = !_is_num(b,i)

@propagate_inbounds function _store_num!(b::RelaxCache{V,N,T}, v::Float64, i::Int) where {V,N,T<:RelaxTag}
    b._num[i] = v
    return
end
@propagate_inbounds function _store_set!(b::RelaxCache{V,N,T}, v::MC{N,T}, i::Int) where {V,N,T<:RelaxTag}
    b._set[i] = v
    return
end
@propagate_inbounds function _store_info!(b::RelaxCache{V,N,T}, v::MC{N,T}, i::Int) where {V,N,T<:RelaxTag}
    b._info[i] = v
    return
end
@propagate_inbounds function _store_subexpression!(b::RelaxCache{V,N,T}, v::MC{N,T}, i::Int) where {V,N,T<:RelaxTag}
    #@inbounds b._subexpression_value[i] = v
    b._subexpression_value[i] = v
    return
end

@inline _first_eval(b::RelaxCache) = b.first_eval
@propagate_inbounds _val(b::RelaxCache{V,N,T}, i::Int) where {V,N,T<:RelaxTag} = _val(b.v, i)
@propagate_inbounds _lbd(b::RelaxCache{V,N,T}, i::Int) where {V,N,T<:RelaxTag} = _lbd(b.v, i)
@propagate_inbounds _ubd(b::RelaxCache{V,N,T}, i::Int) where {V,N,T<:RelaxTag} = _ubd(b.v, i)

_set_input(b::RelaxCache{V,N,T}, n::Int) where {V,N,T<:RelaxTag} = view(b._set_mv_buffer, 1:n)


include(joinpath(@__DIR__, "forward_propagation.jl"))
include(joinpath(@__DIR__, "reverse_propagation.jl"))

function fprop!(t::RELAX_ATTRIBUTE, g::DAT, b::RelaxCache{V,N,T}) where {V,N,T<:RelaxTag}
    f_init!(t, g, b)
    for k = _node_count(g):-1:1
        if _is_unlocked(b, k)
            nt = _node_class(g, k)
            if nt === EXPRESSION
                fprop!(t, Expression(), g, b, k)
            elseif nt === VARIABLE
                fprop!(t, Variable(), g, b, k)
            #elseif nt === SUBEXPRESSION
           #     fprop!(t, Subexpression(), g, b, k)
            end
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
            end
        end
    end
    return flag
end