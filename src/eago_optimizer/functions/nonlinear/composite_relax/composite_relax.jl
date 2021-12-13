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


const RELAX_ATTRIBUTE = Union{Relax,RelaxAA,RelaxAAInfo,RelaxMulEnum,RelaxMulEnumInner,RelaxInterval}
const RELAX_ONLY_ATTRIBUTE = Union{Relax,RelaxAA,RelaxAAInfo,RelaxMulEnum,RelaxMulEnumInner}

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
    ic::IntervalCache                      = IntervalCache{Float64}()
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

    n = node_count(g)
    m = dep_subexpr_count(g)
    p = length(sparsity(g, 1))

    c._num                 = zeros(n)
    c._is_num              = zeros(Bool, n)
    c._set                 = zeros(MC{N,T}, n)
    c._cv_grad_buffer      = zeros(p)
    c._cc_grad_buffer      = zeros(p)
    c._set_mv_buffer       = zeros(MC{N,T}, p)
    c._num_mv_buffer       = zeros(p)
    c._mult_temp           = zero(V)


    for i = 1:n
        ni = node(g, i)
        if node_class(ni) == CONSTANT
            c._is_num[i] = true
            c._num[i] = g.constant_values[first_index(ni)]
        end
    end

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
set(b::RelaxCache{V,N,T}, i::Int) where {V,N,T<:RelaxTag}  = b._set[i]
num(b::RelaxCache{V,N,T}, i::Int) where {V,N,T<:RelaxTag}  = b._num[i]
set_or_num(b::RelaxCache{V,N,T}, i::Int) where {V,N,T<:RelaxTag} = !is_num(b, i) ? set(b,i) : num(b,i) 
info(b::RelaxCache{V,N,T}, i::Int) where {V,N,T<:RelaxTag} = b._info[i]
is_num(b::RelaxCache{V,N,T}) where {V,N,T<:RelaxTag}       = b._is_num
is_num(b::RelaxCache{V,N,T}, i::Int) where {V,N,T<:RelaxTag} = b._is_num[i]
interval(b::RelaxCache{V,N,T}, i::Int) where {V,N,T<:RelaxTag} = Interval{Float64}(set(b, i))
subexpression_set(b::RelaxCache{V,N,T}, i::Int) where {V,N,T<:RelaxTag} = b._subexpression_set[i]
subexpression_num(b::RelaxCache{V,N,T}, i::Int) where {V,N,T<:RelaxTag} = b._subexpression_num[i]
subexpression_is_num(b::RelaxCache{V,N,T}, i::Int) where {V,N,T<:RelaxTag} = b._subexpression_is_num[i]
function Base.setindex!(b::RelaxCache{V,N,T}, v::Float64, i::Int) where {V,N,T<:RelaxTag}
    b._is_num[i] = true
    b._num[i] = v
    nothing
end
function Base.setindex!(b::RelaxCache{V,N,T}, v::MC{N,T}, i::Int) where {V,N,T<:RelaxTag}
    b._is_num[i] = false
    b._set[i] = v
    nothing
end
store_info!(b::RelaxCache{V,N,T}, v::V, i::Int) where {V,N,T<:RelaxTag} = (b._info[i] = v; nothing)

function store_subexpression_set!(b::RelaxCache{V,N,T}, v::MC{N,T}, i::Int) where {V,N,T<:RelaxTag}
    b._subexpression_is_num[i] = false
    b._subexpression_set[i] = v
    nothing
end
function store_subexpression_num!(b::RelaxCache{V,N,T}, v::MC{N,T}, i::Int) where {V,N,T<:RelaxTag}
    b._subexpression_is_num[i] = true
    b._subexpression_num[i] = v
    nothing
end

first_eval(t::RELAX_ONLY_ATTRIBUTE, b::RelaxCache) = b.first_eval
first_eval(t::RelaxInterval, b::RelaxCache) = b.first_eval
val(b::RelaxCache{V,N,T}, i::Int) where {V,N,T<:RelaxTag} = val(b.ic, i)
lbd(b::RelaxCache{V,N,T}, i::Int) where {V,N,T<:RelaxTag} = lbd(b.ic, i)
ubd(b::RelaxCache{V,N,T}, i::Int) where {V,N,T<:RelaxTag} = ubd(b.ic, i)
_set_input(b::RelaxCache{V,N,T}, n::Int) where {V,N,T<:RelaxTag} = view(b._set_mv_buffer, 1:n)
_num_input(b::RelaxCache{V,N,T}, n::Int) where {V,N,T<:RelaxTag} = view(b._num_mv_buffer, 1:n)

include(joinpath(@__DIR__, "utilities.jl"))
include(joinpath(@__DIR__, "forward_propagation.jl"))
include(joinpath(@__DIR__, "reverse_propagation.jl"))

function fprop!(t::RELAX_ATTRIBUTE, g::DAT, b::RelaxCache{V,N,T}) where {V,N,T<:RelaxTag}
    for k = node_count(g):-1:1
        c = node_class(g, k)
        (c == EXPRESSION)    && (fprop!(t, Expression(), g, b, k);    continue)
        (c == VARIABLE)      && (fprop!(t, Variable(), g, b, k);      continue)
        (c == SUBEXPRESSION) && (fprop!(t, Subexpression(), g, b, k); continue)
    end
    nothing
end

function rprop!(t::RELAX_ATTRIBUTE, g::DAT, b::RelaxCache{V,N,T}) where {V,N,T<:RelaxTag}
    flag = r_init!(t, g, b)
    for k = 1:node_count(g)
        nt = node_class(g, k)
        if nt === EXPRESSION
            flag = rprop!(t, Expression(), g, b, k)
        elseif nt === VARIABLE
            flag = rprop!(t, Variable(), g, b, k)
        elseif nt === SUBEXPRESSION
            flag = rprop!(t, Subexpression(), g, b, k)
        end
    end
    return flag
end


function display_expr(g, i)
    n = node(g, i)
    nc = node_class(n)
    if nc == VARIABLE
        return "x[$(first_index(n))]"
    elseif nc == PARAMETER
        return "p[$(first_index(n))]"
    elseif nc == EXPRESSION
        ex_typ = ex_type(n)
        ex_sym = ALL_ATOM_DICT[ex_typ]
        ctup = tuple(children(n)...)
        return "$(ex_sym)$(ctup)"
    elseif nc == CONSTANT
        return "c[$(first_index(n))]"
    end
end

function display_table!(g::DAT, b::RelaxCache{V,N,T}) where {V,N,T<:RelaxTag}
    nc = node_count(g)
    val = [b._is_num[i] ? b._num[i] : b._set[i] for i in 1:nc]
    exr = [display_expr(g, i)  for i in 1:nc]
    data = hcat(exr, b._is_num, val)
    show(pretty_table(data, header = ["Expr", "Is Num", "Val"]; show_row_number = true))
end