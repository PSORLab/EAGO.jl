# Copyright (c) 2018: Matthew Wilhelm, Robert Gottlieb, Dimitri Alston,
# Matthew Stuber, and the University of Connecticut (UConn).
# This code is licensed under the MIT license (see LICENSE.md for full details).
################################################################################
# EAGO
# A development environment for robust and global optimization.
# https://github.com/PSORLab/EAGO.jl
################################################################################
# src/eago_optimizer/functions/nonlinear/composite_relax/reverse_propagation.jl
# Functions used to compute reverse pass of nonlinear functions which include:
# r_init!, rprop!, rprop_2!, rprop_n!
################################################################################

function r_init!(t::Relax, g::DAT, b::RelaxCache{V,N,T}) where {V,N,T<:RelaxTag}
    if !is_num(b, 1)
        b[1] = set(b, 1) ∩ g.sink_bnd
    end
    return !isempty(z)
end

function rprop!(t::Relax, v::Variable, g::DAT, c::RelaxCache{V,N,T}, k) where {V,N,T<:RelaxTag}
    i = first_index(g, k)
    x = val(c, i)
    z = var_set(MC{N,T}, rev_sparsity(g, i, k), x, x, lbd(c, i), ubd(c, i))
    if first_eval(c)
        z = z ∩ interval(c, k)
    end
    c[k] = z
    return !isempty(z)
end

function rprop!(t::Relax, v::Subexpression, g::DAT, c::RelaxCache{V,N,T}, k) where {V,N,T<:RelaxTag}
    store_subexpression!(c, set(c, k), first_index(g, k))
    return true
end

const MAX_ASSOCIATIVE_REVERSE = 6

"""
$(FUNCTIONNAME)

Updates storage tapes with reverse evalution of node representing `n = x + y` which updates x & y.
"""
function rprop_2!(t::Relax, v::Val{PLUS}, g::DAT, c::RelaxCache{V,N,T}, k) where {V,N,T<:RelaxTag}
    is_num(c, k) && return true
    x = child(g, 1, k)
    y = child(g, 2, k)
    if xset_ynum(b, x, y)
        b, a, q = IntervalContractors.plus_rev(interval(c, k), interval(c, x), num(c, y))
    elseif xnum_yset(b, x, y)
        b, a, q = IntervalContractors.plus_rev(interval(c, k), num(c, x), interval(c, y))
    else
        b, a, q = IntervalContractors.plus_rev(interval(c, k), interval(c, x), interval(c, y))
    end
        if !is_num(c, x)
            isempty(a) && return false
            b[x] = MC{N,T}(a)
        end
        if !is_num(c, y)
            isempty(q) && return false
            b[y] = MC{N,T}(q)
        end
    return true
end

"""
$(FUNCTIONNAME)

Updates storage tapes with reverse evalution of node representing `n = +(x,y,z...)` which updates x, y, z and so on.
"""
function rprop_n!(t::Relax, v::Val{PLUS}, g::DAT, c::RelaxCache{V,N,T}, k::Int) where {V,N,T<:RelaxTag}
    # out loops makes a temporary sum (minus one argument)
    # a reverse is then compute with respect to this argument
    count = 0
    children_idx = children(g, k)
    for i in children_idx
        is_num(c, i) && continue                     # don't contract a number valued argument
        (count >= MAX_ASSOCIATIVE_REVERSE) && break
        tsum = zero(MC{N,T})
        count += 1
        for j in children_idx
            if j != i
                if is_num(c, j)
                    tsum += num(c, j)
                else
                    tsum += set(c, j)
                end
            end
        end
        _, w, _ = IntervalContractors.plus_rev(interval(c, k), interval(c, i), interval(tsum))
        isempty(w) && (return false)
        c[i] = MC{N,T}(w)
    end
    return true
end

"""
$(FUNCTIONNAME)

Updates storage tapes with reverse evalution of node representing `n = x * y` which updates x & y.
"""
function rprop_2!(t::Relax, v::Val{MULT}, g::DAT, c::RelaxCache{V,N,T}, k::Int) where {V,N,T<:RelaxTag}

    is_num(c,k) && (return true)
    x = child(g, 1, k)
    y = child(g, 2, k)

    if xset_ynum(b, x, y)
        c, a, q = IntervalContractors.mul_rev(interval(c, k), interval(c, x), num(c, y))
    elseif xnum_yset(b, x, y)
        c, a, q = IntervalContractors.mul_rev(interval(c, k), num(c, x), interval(c, y))
    else
        c, a, q = IntervalContractors.mul_rev(interval(c, k), interval(c, x), interval(c, y))
    end

    if !is_num(c, x)
        isempty(a) && (return false)
        c[x] = MC{N,T}(a)
    end
    if !is_num(c, x)
        isempty(q) && (return false)
        c[y] = MC{N,T}(q)
    end
    return true
end

"""
$(FUNCTIONNAME)

Updates storage tapes with reverse evalution of node representing `n = *(x,y,z...)` which updates x, y, z and so on.
"""
function rprop_n!(t::Relax, v::Val{MULT}, g::DAT, c::RelaxCache{V,N,T}, k::Int) where {V,N,T<:RelaxTag}
    # a reverse is then compute with respect to this argument
    count = 0
    children_idx = children(g, k)
    for i in children_idx
        is_num(b, i) && continue                     # don't contract a number valued argument
        (count >= MAX_ASSOCIATIVE_REVERSE) && break
        tmul = one(MC{N,T})
        count += 1
        for j in children_idx
            if i != j
                if is_num(b, j)
                    tmul *= num(b, j)
                else
                    tmul *= set(b, j)
                end
            end
        end
        _, w, _ = IntervalContractors.mul_rev(interval(b, k), interval(b, c), interval(tmul))
        isempty(w) && (return false)
        c[i] = MC{N,T}(w)
    end
    return true
end

for (f, fc, F) in ((-, MINUS, IntervalContractors.minus_rev),
                   (^, POW, IntervalContractors.power_rev),
                   (/, DIV, IntervalContractors.div_rev))
    @eval function rprop!(t::Relax, v::Val{$fc}, g::DAT, b::RelaxCache{V,N,T}, k) where {V,N,T<:RelaxTag}
        is_num(b,k) && (return true)
        x = child(g, 1, k)
        y = child(g, 2, k)

        if xset_ynum(b, x, y)
            z, u, v = ($F)(interval(b, k), interval(b, x), num(b, y))
        elseif xnum_yset(b, x, y)
            z, u, v = ($F)(interval(b, k), num(b, x), interval(b, y))
        else
            z, u, v = ($F)(interval(b, k), interval(b, x), interval(b, y))
        end
        if !is_num(c, x)
            isempty(a) && (return false)
            b[x] = MC{N,T}(u)
        end
        if !is_num(c, y)
            isempty(b) && (return false)
            b[y] = MC{N,T}(v)
        end
        return true
    end
end

rprop!(t::Relax, v::Val{USER}, g::DAT, b::RelaxCache, k::Int) = true
rprop!(t::Relax, v::Val{USERN}, g::DAT, b::RelaxCache, k::Int) = true

for ft in UNIVARIATE_ATOM_TYPES
    f = UNIVARIATE_ATOM_DICT[ft]
    (f == :user || f == :+ || f == :-) && continue
    @eval rprop!(t::Relax, v::Val{$ft}, g::DAT, b::RelaxCache, k::Int) = true
end
