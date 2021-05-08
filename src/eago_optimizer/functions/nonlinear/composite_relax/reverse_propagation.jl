# Copyright (c) 2018: Matthew Wilhelm & Matthew Stuber.
# This code is licensed under MIT license (see LICENSE.md for full details)
#############################################################################
# EAGO
# A development environment for robust and global optimization
# See https://github.com/PSORLab/EAGO.jl
#############################################################################
# src/eago_optimizer/functions/nonlinear/reverse_pass.jl
# Functions used to compute reverse pass of nonlinear functions.
#############################################################################

function rprop!(::Type{Relax}, ::Type{Variable}, g::AbstractDG, b::RelaxCache{V,S}, k::Int) where {V,S}
    i = _first_index(g, k)
    x = _val(b, i)
    z = _var_set(V, _rev_sparsity(g, i), x, x, _lbd(b, i), _ubd(b, i))
    if _first_eval(b)
        z = z ∩ _interval(b, k)
    end
    _store_set!(b, z, k)
    return
end

function rprop!(::Type{Relax}, ::Type{Subexpression}, g::AbstractDG, c::RelaxCache{V,S}, k::Int) where {V,S}
    _store_subexpression!(c, _set(c, k), _first_index(g, k))
end

const MAX_ASSOCIATIVE_REVERSE = 6

"""
$(FUNCTIONNAME)

Updates storage tapes with reverse evalution of node representing `n = x + y` which updates x & y.
"""
function rprop_2!(::Type{Relax}, ::typeof(+), g::AbstractDG, b::RelaxCache{V,S}, k::Int) where {V,S}

    _is_num(b, k) && (return true)
    x = _child(g, 1, k)
    y = _child(g, 2, k)
    x_is_num = _is_num(b, x)
    y_is_num = _is_num(b, x)

    if !x_is_num && y_is_num
        c, a, q = IntervalContractors.plus_rev(_interval(b, k), _interval(b, x), _num(b, y))
    elseif x_is_num && !y_is_num
        c, a, q = IntervalContractors.plus_rev(_interval(b, k), _num(b, x), _interval(b, y))
    else
        c, a, q = IntervalContractors.plus_rev(_interval(b, k), _interval(b, x), _interval(b, y))
    end

    if !x_is_num
        isempty(a) && (return false)
        _store_set!(b, V(a), x)
    end
    if !y_is_num
        isempty(q) && (return false)
        _store_set!(b, V(q), y)
    end
    return true
end

"""
$(FUNCTIONNAME)

Updates storage tapes with reverse evalution of node representing `n = +(x,y,z...)` which updates x, y, z and so on.
"""
function rprop_n!(::Type{Relax}, ::typeof(+), g::AbstractDG, b::RelaxCache{V,S}, k::Int) where {V,S}
    # out loops makes a temporary sum (minus one argument)
    # a reverse is then compute with respect to this argument
    count = 0
    children_idx = _children(g, k)
    for c in children_idx
        is_num(b, c) && continue                     # don't contract a number valued argument
        (count >= MAX_ASSOCIATIVE_REVERSE) && break
        tsum = zero(V)
        count += 1
        for i in children_idx
            if i != c
                if numvalued[i]
                    tsum += _num(b, i)
                else
                    tsum += _set(b, i)
                end
            end
        end
        q, w, z = IntervalContractors.plus_rev(_interval(b, k), _interval(b, c), interval(tsum))
        isempty(w) && (return false)
        _store_set!(b, V(w), c)
    end
    return true
end

"""
$(FUNCTIONNAME)

Updates storage tapes with reverse evalution of node representing `n = x * y` which updates x & y.
"""
function rprop_2!(::Type{Relax}, ::typeof(*), g::AbstractDG, b::RelaxCache{V,S}, k::Int) where {V,S}

    _is_num(b,k) && (return true)
    x = _child(g, 1, k)
    y = _child(g, 2, k)
    x_is_num = _is_num(b, x)
    y_is_num = _is_num(b, x)

    if !x_is_num && y_is_num
        c, a, q = IntervalContractors.mul_rev(_interval(b, k), _interval(b, x), _num(b, y))
    elseif x_is_num && !y_is_num
        c, a, q = IntervalContractors.mul_rev(_interval(b, k), _num(b, x), _interval(b, y))
    else
        c, a, q = IntervalContractors.mul_rev(_interval(b, k), _interval(b, x), _interval(b, y))
    end

    if !x_is_num
        isempty(a) && (return false)
        _store_set!(b, V(a), x)
    end
    if !y_is_num
        isempty(q) && (return false)
        _store_set!(b, V(q), y)
    end
    return true
end

"""
$(FUNCTIONNAME)

Updates storage tapes with reverse evalution of node representing `n = *(x,y,z...)` which updates x, y, z and so on.
"""
function rprop_n!(::Type{Relax}, ::typeof(*), g::AbstractDG, b::RelaxCache{V,S}, k::Int) where {V,S}
    # a reverse is then compute with respect to this argument
    count = 0
    children_idx = _children(g, k)
    for c in children_idx
        is_num(b, c) && continue                     # don't contract a number valued argument
        (count >= MAX_ASSOCIATIVE_REVERSE) && break
        tmul = one(V)
        count += 1
        for i in children_idx
            if i != c
                if _is_num(b, i)
                    tmul *= _num(b, i)
                else
                    tmul *= _set(b, i)
                end
            end
        end
        q, w, z = IntervalContractors.mul_rev(_interval(b, k), _interval(b, c), interval(tmul))
        isempty(w) && (return false)
        _store_set!(b, V(w), c)
    end
    return true
end

for (f, F) in ((-, IntervalContractors.minus_rev), (^, IntervalContractors.power_rev), (/, IntervalContractors.div_rev))
    @eval function rprop!(::Type{Relax}, ::typeof($f), g::AbstractDG, b::RelaxCache{V,S}, k::Int) where {V,S}
        _is_num(b,k) && (return true)
        x = _child(g, 1, k)
        y = _child(g, 2, k)
        x_is_num = _is_num(b, x)
        y_is_num = _is_num(b, x)

        if !x_is_num && y_is_num
            z, u, v = ($F)(_interval(b, k), _interval(b, x), _num(b, y))
        elseif x_is_num && y_is_num
            z, u, v = ($F)(_interval(b, k), _num(b, x), _interval(b, y))
        else
            z, u, v = ($F)(_interval(b, k), _interval(b, x), _interval(b, y))
        end

        if !x_is_num
            isempty(a) && (return false)
            _store_set!(b, V(u), x)
        end
        if !y_is_num
            isempty(b) && (return false)
            _store_set!(b, V(v), y)
        end
        return true
    end
end

rprop!(::Type{Relax}, ::typeof(user), g::AbstractDG, b::RelaxCache, k::Int) = nothing
rprop!(::Type{Relax}, ::typeof(usern), g::AbstractDG, b::RelaxCache, k::Int) = nothing

# TODO: Define individual reverse univariates...
