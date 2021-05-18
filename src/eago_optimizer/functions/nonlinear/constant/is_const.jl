#=
Routines used to prepopulate numeric expressions and parameters in Set-valued buffers
We probably don't need a separate buffer for each function rather we should just
use a single one of the largest size expression.
=#
struct ConstInfo <: AbstractCacheAttribute end
struct ConstData{T<:Real}
    val::T
    is_num::Bool
end
ConstData(val::T, b::Bool) where T<:Real = ConstData{T}(x, b)
ConstData(::T) where T<:Real = ConstData(zero(T), false)
@inline _val(x::ConstData{T}) where {T<:Real} = x.val
@inline _is_num(x::ConstData) = x.is_num

mutable struct ConstCache{T} <: AbstractCache
    val::Vector{T}
    is_num::Vector{Bool}
end
ConstCache(::T, n::Int) where T<:Real = ConstCache{T}(fill(zero(T), n), fill(false, n))
@inline _val(b::ConstCache{T}, i::Int) where T = getindex(b.val, i)
@inline _is_num(b::ConstCache{T}, i::Int) where T = getindex(b.is_num, i)
@inline function _store_data!(b::ConstCache{T}, i::Int, v::ConstData{T}) where T
    b.val[i] = v.val
    b.is_num[i] = v.is_num
    return
end

for (F,C) in ((:+, PLUS), (:*, MULT), (:min, MIN), (:max, MAX), (:usern, USERN))
    @eval function fprop!(::ConstInfo, ::Val{$F}, g::DAG, b::ConstCache{T}, k::Int) where T
        clist = _children(g, k)
        if any(i -> !_is_num(b, i), clist)
            z = ConstData(T)
        else
            z = ConstData(mapreduce(i -> _val(b, i), $C, clist), false)
        end
        _store_data!(b, k, z)
    end
end

function fprop!(::Type{ConstInfo}, ::Val{POW}, g::DAG, b::ConstCache{T}, k::Int) where T
    x = _child(g, 1, k)
    y = _child(g, 2, k)
    xc = _val(b, x)
    yc = _val(b, y)
    if (_is_num(b, x) && _is_num(b, y)) || iszero(yc)
        z = ConstData(xc^yc, true)
    else
        z = ConstData(T)
    end
    store_data!(b, k, z)
end

for F in (/, -, arh, lower_bnd, upper_bnd)
    @eval function fprop!(::Type{ConstInfo}, ::typeof($F), g::DAG, b::ConstCache{T}, k::Int) where T
        xc = _data(b, _child(g, 1, k))
        yc = _data(b, _child(g, 2, k))
        if _is_num(xc) && _is_num(yc)
            z = ConstData(($F)(_val(xc),_val(yc)), true)
        else
            z = ConstData(T)
        end
        store_data!(b, k, z)
    end
end

for F in values(UNIVARIATE_EVAL)
    @eval function fprop!(::Type{ConstInfo}, ::typeof($F), g::DAG, b::ConstCache{T}, k::Int) where T
        xc = _data(b, _first_child(g, k))
        z = _is_num(xc) ? ConstData(T) : ConstData(($F)(_val(xc)), true)
        store_data!(b, k, z)
    end
end
