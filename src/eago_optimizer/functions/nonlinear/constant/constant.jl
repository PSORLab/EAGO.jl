#=
Routines used to prepopulate numeric expressions and parameters in Set-valued buffers
We probably don't need a separate buffer for each function rather we should just
use a single one of the largest size expression.
=#
struct ConstInfo <: AbstractCacheAttribute end
struct ConstData
    _num::Float64
    _is_num::Bool
end
ConstData() = ConstData(0.0, false)
_num(x::ConstData) = x._num
_is_num(x::ConstData) = x._is_num

mutable struct ConstantCache <: AbstractCache
    _num::Vector{Float64}
    _is_num::Vector{Bool}
end
ConstantCache(n::Int) = ConstantCache(zeros(n), fill(false, n))

_num(b::ConstantCache) = b._num
_num(b::ConstantCache, i::Int) = getindex(b._num, i)
_is_num(b::ConstantCache) = b._is_num
_is_num(b::ConstantCache, i::Int) = getindex(b._is_num, i)

function _store_data!(b::ConstantCache, i::Int, v::ConstData)
    b._num[i] = v._num
    b._is_num[i] = v._is_num
    return
end

function initialize!(c::ConstantCache, g::ALLGRAPHS)
    n = _node_count(g)
    c._num     = zeros(Float64, n)
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
eval(quote
    function fprop!(::ConstInfo, ::Val{$C}, g::ALLGRAPHS, b::ConstantCache, k::Int)
        clist = _children(g, k)
        if any(i -> !_is_num(b, i), clist)
            z = ConstData()
        else
            z = ConstData(mapreduce(i -> _num(b, i), $F, clist), false)
        end
        _store_data!(b, k, z)
    end
end)
end

function fprop_2!(::ConstInfo, ::Val{MINUS}, g::ALLGRAPHS, b::ConstantCache, k::Int)
    x = _child(g, 1, k)
    y = _child(g, 2, k)
    if _is_num(b, x) && _is_num(b, y)
        z = ConstData(_num(b, x) - _num(b, y), true)
    else
        z = ConstData()
    end
    store_data!(b, k, z)
end
function fprop_n!(::ConstInfo, ::Val{MINUS}, g::ALLGRAPHS, b::ConstantCache, k::Int)
    x = _child(g, 1, k)
    y = _child(g, 2, k)
    if _is_num(b, x) && _is_num(b, y)
        z = ConstData(_num(b, x) - _num(b, y), true)
    else
        z = ConstData()
    end
    store_data!(b, k, z)
end
function fprop!(::ConstInfo, ::Val{MINUS}, g::ALLGRAPHS, b::ConstantCache, k::Int)
    n = _arity(g, k)
    if n == 2
        return fprop_2!(ConstInfo(), Val(MINUS), g, b, k)
    end
    fprop_n!(ConstInfo(), Val(MINUS), g, b, k)
end

function fprop!(::ConstInfo, ::Val{POW}, g::ALLGRAPHS, b::ConstantCache, k::Int)
    x = _child(g, 1, k)
    y = _child(g, 2, k)
    xc = _num(b, x)
    yc = _num(b, y)
    if (_is_num(b, x) && _is_num(b, y)) || iszero(yc)
        z = ConstData(xc^yc, true)
    else
        z = ConstData()
    end
    store_data!(b, k, z)
end

for F in BIVARIATE_ATOM_TYPES
    f = BIVARIATE_ATOM_DICT[F]
    if f == :- || f == :+ || f == :* || f == :min || f == :max || f == :usern || f == :^
        continue
    end
    eval(quote 
        function fprop!(::ConstInfo, ::Val{$F}, g::ALLGRAPHS, b::ConstantCache, k::Int)
            x = _child(g, 1, k)
            y = _child(g, 2, k)
            if _is_num(b, x) && _is_num(b, y)
                z = ConstData(($f)(_num(b, x),_num(b, y)), true)
            else
                z = ConstData(T)
            end
            store_data!(b, k, z)
        end
    end)
end

for F in UNIVARIATE_ATOM_TYPES
    f = UNIVARIATE_ATOM_DICT[F]
    if f == :- || f == :+ || f == :* || f == :min || f == :max || f == :usern
        continue
    end
    eval(quote 
        function fprop!(::ConstInfo, ::Val{$F}, g::ALLGRAPHS, b::ConstantCache, k::Int)
            x = _first_child(g, k)
            z = _is_num(b, x) ? ConstData() : ConstData(($f)(_num(b, x)), true)
            store_data!(b, k, z)
        end
    end)
end


function fprop!(::ConstInfo, g::ALLGRAPHS, b::ConstantCache)
    f_init!(t, g, b)
    for k = _node_count(g):-1:1
        if _is_unlocked(b, k)
            nt = _node_class(g, k)
            if nt == EXPRESSION
                fprop!(ConstInfo(), Expression(), g, b, k)
            elseif nt == VARIABLE
                fprop!(ConstInfo(), Variable(), g, b, k)
            elseif nt == SUBEXPRESSION
                fprop!(ConstInfo(), Subexpression(), g, b, k)
            end
        end
    end
    return
end

function rprop!(::ConstInfo, g::ALLGRAPHS, b::ConstantCache)
    flag = r_init!(t, g, b)
    for k = 1:_node_count(g)
        if _is_unlocked(b, k)
            nt = _node_class(g, k)
            if nt == EXPRESSION
                flag = rprop!(ConstInfo(), Expression(), g, b, k)
            elseif nt == VARIABLE
                flag = rprop!(ConstInfo(), Variable(), g, b, k)
            end
        end
    end
    return flag
end