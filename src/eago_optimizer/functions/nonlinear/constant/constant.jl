#=
Routines used to prepopulate numeric expressions and parameters in Set-valued buffers
We probably don't need a separate buffer for each function rather we should just
use a single one of the largest size expression.
=#
struct ConstInfo <: AbstractCacheAttribute end
struct ConstData{T<:Real}
    _num::T
    _is_num::Bool
end
ConstData(::T) where T<:Real = ConstData(zero(T), false)
@inline _num(x::ConstData{T}) where {T<:Real} = x._num
@inline _is_num(x::ConstData) = x._is_num

mutable struct ConstantCache{T} <: AbstractCache
    _num::Vector{T}
    _is_num::Vector{Bool}
end
ConstantCache(::T, n::Int) where T<:Real = ConstantCache{T}(fill(zero(T), n), fill(false, n))
ConstantCache{T}(n::Int) where T<:Real = ConstantCache(zero(T), n)

@inline _num(b::ConstantCache{T}) where {T<:Real} = b._num
@inline _num(b::ConstantCache{T}, i::Int) where {T<:Real} = getindex(b._num, i)
@inline _is_num(b::ConstantCache) = b._is_num
@inline _is_num(b::ConstantCache, i::Int) = getindex(b._is_num, i)

@inline function _store_data!(b::ConstantCache{T}, i::Int, v::ConstData{T}) where T
    b._num[i] = v._num
    b._is_num[i] = v._is_num
    return
end

function initialize!(c::ConstantCache{S}, g::AbstractDG) where S<:Real
    n = _node_count(g)
    c._num     = zeros(S, n)
    c._is_num  = zeros(Bool, n)
    for i = 1:n
        nd = _node(g, i)
        if _ex_type(nd) == CONST_ATOM
            c._is_num[i] = true
            c._num[i] = _constant_value(g, _first_index(nd))
        elseif _ex_type(nd) == PARAM_ATOM
            c._is_num[i] = true
            c._num[i] = _parameter_value(g, _first_index(nd))
        end
    end
    return
end

for (F,C) in ((:+, PLUS), (:*, MULT), (:min, MIN), (:max, MAX), (:usern, USERN))
    @eval function fprop!(::ConstInfo, ::Val{$C}, g::AbstractDG, b::ConstantCache{T}, k::Int) where T
        clist = _children(g, k)
        if any(i -> !_is_num(b, i), clist)
            z = ConstData(T)
        else
            z = ConstData(mapreduce(i -> _num(b, i), $F, clist), false)
        end
        _store_data!(b, k, z)
    end
end

function fprop_2!(::ConstInfo, ::Val{MINUS}, g::AbstractDG, b::ConstantCache{T}, k::Int) where T
    x = _child(g, 1, k)
    y = _child(g, 2, k)
    if _is_num(b, x) && _is_num(b, y)
        z = ConstData(_num(b, x) - _num(b, y), true)
    else
        z = ConstData(T)
    end
    store_data!(b, k, z)
end
function fprop_n!(::ConstInfo, ::Val{MINUS}, g::AbstractDG, b::ConstantCache{T}, k::Int) where T
    x = _child(g, 1, k)
    y = _child(g, 2, k)
    if _is_num(b, x) && _is_num(b, y)
        z = ConstData(_num(b, x) - _num(b, y), true)
    else
        z = ConstData(T)
    end
    store_data!(b, k, z)
end
function fprop!(::ConstInfo, ::Val{MINUS}, g::AbstractDG, b::ConstantCache{T}, k::Int) where T
    n = _arity(g, k)
    if n == 2
        return fprop_2!(ConstInfo(), Val(MINUS), g, b, k)
    end
    fprop_n!(ConstInfo(), Val(MINUS), g, b, k)
end

function fprop!(::ConstInfo, ::Val{POW}, g::AbstractDG, b::ConstantCache{T}, k::Int) where T
    x = _child(g, 1, k)
    y = _child(g, 2, k)
    xc = _num(b, x)
    yc = _num(b, y)
    if (_is_num(b, x) && _is_num(b, y)) || iszero(yc)
        z = ConstData(xc^yc, true)
    else
        z = ConstData(T)
    end
    store_data!(b, k, z)
end

for F in BIVARIATE_ATOM_TYPES
    f = BIVARIATE_ATOM_DICT[F]
    if f == :- || f == :+ || f == :* || f == :min || f == :max || f == :usern || f == :^
        continue
    end
    @eval function fprop!(::ConstInfo, ::Val{$F}, g::AbstractDG, b::ConstantCache{T}, k::Int) where T
        x = _child(g, 1, k)
        y = _child(g, 2, k)
        if _is_num(b, x) && _is_num(b, y)
            z = ConstData(($f)(_num(b, x),_num(b, y)), true)
        else
            z = ConstData(T)
        end
        store_data!(b, k, z)
    end
end

for F in UNIVARIATE_ATOM_TYPES
    f = UNIVARIATE_ATOM_DICT[F]
    if f == :- || f == :+ || f == :* || f == :min || f == :max || f == :usern
        continue
    end
    @eval function fprop!(::ConstInfo, ::Val{$F}, g::AbstractDG, b::ConstantCache{T}, k::Int) where T
        x = _first_child(g, k)
        z = _is_num(b, x) ? ConstData(T) : ConstData(($f)(_num(b, x)), true)
        store_data!(b, k, z)
    end
end
