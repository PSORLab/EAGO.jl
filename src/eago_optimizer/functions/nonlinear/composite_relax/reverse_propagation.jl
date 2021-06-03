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

function r_init!(t::Relax, g::DAT, c::RelaxCache{N,T}) where {N,T<:RelaxTag}
    if !_is_num(c, 1)
        z = _set(c, 1) ∩ g.sink_bnd
        _store_set!(c, z, 1)
    end
    return !isempty(z)
end

function rprop!(t::Relax, v::Variable, g::DAT, c::RelaxCache{N,T}, k::Int) where {N,T<:RelaxTag}
    i = _first_index(g, k)
    x = _val(c, i)
    z = _var_set(MC{N,T}, _rev_sparsity(g, i, k), x, x, _lbd(c, i), _ubd(c, i))
    if _first_eval(c)
        z = z ∩ _interval(c, k)
    end
    _store_set!(c, z, k)
    return !isempty(z)
end

function rprop!(t::Relax, v::Subexpression, g::DAT, c::RelaxCache{N,T}, k::Int) where {N,T<:RelaxTag}
    _store_subexpression!(c, _set(c, k), _first_index(g, k))
    return true
end

const MAX_ASSOCIATIVE_REVERSE = 6

"""
$(FUNCTIONNAME)

Updates storage tapes with reverse evalution of node representing `n = x + y` which updates x & y.
"""
function rprop_2!(t::Relax, v::Val{PLUS}, g::DAT, c::RelaxCache{N,T}, k::Int) where {N,T<:RelaxTag}

    _is_num(c, k) && (return true)
    x = _child(g, 1, k)
    y = _child(g, 2, k)
    x_is_num = _is_num(c, x)
    y_is_num = _is_num(c, x)

    if !x_is_num && y_is_num
        b, a, q = IntervalContractors.plus_rev(_interval(c, k), _interval(c, x), _num(c, y))
    elseif x_is_num && !y_is_num
        b, a, q = IntervalContractors.plus_rev(_interval(c, k), _num(c, x), _interval(c, y))
    else
        b, a, q = IntervalContractors.plus_rev(_interval(c, k), _interval(c, x), _interval(c, y))
    end

    if !x_is_num
        isempty(a) && (return false)
        _store_set!(b, MC{N,T}(a), x)
    end
    if !y_is_num
        isempty(q) && (return false)
        _store_set!(b, MC{N,T}(q), y)
    end
    return true
end

"""
$(FUNCTIONNAME)

Updates storage tapes with reverse evalution of node representing `n = +(x,y,z...)` which updates x, y, z and so on.
"""
function rprop_n!(t::Relax, v::Val{PLUS}, g::DAT, c::RelaxCache{N,T}, k::Int) where {N,T<:RelaxTag}
    # out loops makes a temporary sum (minus one argument)
    # a reverse is then compute with respect to this argument
    count = 0
    children_idx = _children(g, k)
    for i in children_idx
        _is_num(c, i) && continue                     # don't contract a number valued argument
        (count >= MAX_ASSOCIATIVE_REVERSE) && break
        tsum = zero(MC{N,T})
        count += 1
        for j in children_idx
            if j != i
                if _is_num(c, j)
                    tsum += _num(c, j)
                else
                    tsum += _set(c, j)
                end
            end
        end
        _, w, _ = IntervalContractors.plus_rev(_interval(c, k), _interval(c, i), interval(tsum))
        isempty(w) && (return false)
        _store_set!(c, MC{N,T}(w), i)
    end
    return true
end

"""
$(FUNCTIONNAME)

Updates storage tapes with reverse evalution of node representing `n = x * y` which updates x & y.
"""
function rprop_2!(t::Relax, v::Val{MULT}, g::DAT, c::RelaxCache{N,T}, k::Int) where {N,T<:RelaxTag}

    _is_num(b,k) && (return true)
    x = _child(g, 1, k)
    y = _child(g, 2, k)
    x_is_num = _is_num(c, x)
    y_is_num = _is_num(c, y)

    if !x_is_num && y_is_num
        c, a, q = IntervalContractors.mul_rev(_interval(c, k), _interval(c, x), _num(c, y))
    elseif x_is_num && !y_is_num
        c, a, q = IntervalContractors.mul_rev(_interval(c, k), _num(c, x), _interval(c, y))
    else
        c, a, q = IntervalContractors.mul_rev(_interval(c, k), _interval(c, x), _interval(c, y))
    end

    if !x_is_num
        isempty(a) && (return false)
        _store_set!(c, MC{N,T}(a), x)
    end
    if !y_is_num
        isempty(q) && (return false)
        _store_set!(c, MC{N,T}(q), y)
    end
    return true
end

"""
$(FUNCTIONNAME)

Updates storage tapes with reverse evalution of node representing `n = *(x,y,z...)` which updates x, y, z and so on.
"""
function rprop_n!(t::Relax, v::Val{MULT}, g::DAT, c::RelaxCache{N,T}, k::Int) where {N,T<:RelaxTag}
    # a reverse is then compute with respect to this argument
    count = 0
    children_idx = _children(g, k)
    for i in children_idx
        _is_num(b, i) && continue                     # don't contract a number valued argument
        (count >= MAX_ASSOCIATIVE_REVERSE) && break
        tmul = one(MC{N,T})
        count += 1
        for j in children_idx
            if i != j
                if _is_num(b, j)
                    tmul *= _num(b, j)
                else
                    tmul *= _set(b, j)
                end
            end
        end
        _, w, _ = IntervalContractors.mul_rev(_interval(b, k), _interval(b, c), interval(tmul))
        isempty(w) && (return false)
        _store_set!(c, MC{N,T}(w), i)
    end
    return true
end

for (f, fc, F) in ((-, MINUS, IntervalContractors.minus_rev),
                   (^, POW, IntervalContractors.power_rev),
                   (/, DIV, IntervalContractors.div_rev))
    eval(quote 
        function rprop!(t::Relax, v::Val{$fc}, g::DAT, b::RelaxCache{N,T}, k::Int) where {N,T<:RelaxTag}
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
                _store_set!(b, MC{N,T}(u), x)
            end
            if !y_is_num
                isempty(b) && (return false)
                _store_set!(b, MC{N,T}(v), y)
            end
            return true
        end
    end)
end

rprop!(t::Relax, v::Val{USER}, g::DAT, b::RelaxCache, k::Int) = true
rprop!(t::Relax, v::Val{USERN}, g::DAT, b::RelaxCache, k::Int) = true

for ft in UNIVARIATE_ATOM_TYPES
    f = UNIVARIATE_ATOM_DICT[ft]
    if f == :user || f == :+ || f == :-
        continue
    end
    @eval rprop!(t::Relax, v::Val{$ft}, g::DAT, b::RelaxCache, k::Int) = true
end
