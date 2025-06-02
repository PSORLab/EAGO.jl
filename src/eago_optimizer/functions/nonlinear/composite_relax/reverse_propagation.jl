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
    return !isempty(set(b, 1))
end

function rprop!(t::RelaxCacheAttribute, v::Variable, g::DAT, b::RelaxCache{V,N,T}, k) where {V,N,T<:RelaxTag}
    i = first_index(g, k)
    x = val(b, i)
    l = lbd(b, i)
    u = ubd(b, i)
    if l == u
        b[k] = x
    else
        z = varset(MC{N,T}, rev_sparsity(g, i, k), x, x, l, u)
        if first_eval(t, b)
            z = z ∩ interval(b, k)
        end
        b[k] = z
        b.ic.v.lower_variable_bounds[i] = z.Intv.lo
        b.ic.v.upper_variable_bounds[i] = z.Intv.hi
    end
    return !isempty(set(b, k))
end

function rprop!(t::Relax, ex::Subexpression, g::DAT, b::RelaxCache{V,N,T}, k) where {V,N,T<:RelaxTag}
    i = first_index(g, k)
    x = dependent_subexpression_index(g, i)
    if subexpression_is_num(b, x)
        b[k] = subexpression_num(b, x)
    else
        b[k] = subexpression_set(b, x)
    end
    return true
end

for F in (PLUS, MULT)
    @eval function rprop!(t::Relax, v::Val{$F}, g::DAT, b::RelaxCache{V,N,T}, k) where {V,N,T<:RelaxTag}
        is_binary(g, k) ? rprop_2!(Relax(), Val($F), g, b, k) : rprop_n!(Relax(), Val($F), g, b, k)
    end
end

const MAX_ASSOCIATIVE_REVERSE = 6

"""
$(FUNCTIONNAME)

Updates storage tapes with reverse evaluation of node representing `n = x + y` which updates x and y.
"""
function rprop_2!(t::Relax, v::Val{PLUS}, g::DAT, c::RelaxCache{V,N,T}, k) where {V,N,T<:RelaxTag}
    is_num(c, k) && return true
    x = child(g, 1, k)
    y = child(g, 2, k)
    if xset_ynum(c, x, y)
        b, a, q = IntervalContractors.plus_rev(interval(c, k), interval(c, x), num(c, y))
    elseif xnum_yset(c, x, y)
        b, a, q = IntervalContractors.plus_rev(interval(c, k), num(c, x), interval(c, y))
    else
        b, a, q = IntervalContractors.plus_rev(interval(c, k), interval(c, x), interval(c, y))
    end
        if !is_num(c, x)
            isempty(a) && return false
            c[x] = MC{N,T}(a)
        end
        if !is_num(c, y)
            isempty(q) && return false
            c[y] = MC{N,T}(q)
        end
    return true
end

"""
$(FUNCTIONNAME)

Updates storage tapes with reverse evaluation of node representing `n = +(x,y,z...)` which updates x, y, z and so on.
"""
function rprop_n!(t::Relax, v::Val{PLUS}, g::DAT, c::RelaxCache{V,N,T}, k::Int) where {V,N,T<:RelaxTag}
    # Outer loop makes a temporary sum (minus one argument)
    # A reverse is then computed with respect to this argument
    count = 0
    children_idx = children(g, k)
    for i in children_idx
        is_num(c, i) && continue                     # Don't contract a number valued argument
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
        _, w, _ = IntervalContractors.plus_rev(interval(c, k), interval(c, i), Intv(tsum))
        isempty(w) && (return false)
        c[i] = MC{N,T}(w)
    end
    return true
end

"""
$(FUNCTIONNAME)

Updates storage tapes with reverse evaluation of node representing `n = x * y` which updates x and y.
"""
function rprop_2!(t::Relax, v::Val{MULT}, g::DAT, c::RelaxCache{V,N,T}, k::Int) where {V,N,T<:RelaxTag}

    is_num(c, k) && (return true)
    x = child(g, 1, k)
    y = child(g, 2, k)

    if xset_ynum(c, x, y)
        b, a, q = IntervalContractors.mul_rev(interval(c, k), interval(c, x), num(c, y))
    elseif xnum_yset(c, x, y)
        b, a, q = IntervalContractors.mul_rev(interval(c, k), num(c, x), interval(c, y))
    else
        b, a, q = IntervalContractors.mul_rev(interval(c, k), interval(c, x), interval(c, y))
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

Updates storage tapes with reverse evaluation of node representing `n = *(x,y,z...)` which updates x, y, z and so on.
"""
function rprop_n!(t::Relax, v::Val{MULT}, g::DAT, b::RelaxCache{V,N,T}, k::Int) where {V,N,T<:RelaxTag}
    # Outer loop makes a temporary sum (minus one argument)
    # A reverse is then computed with respect to this argument
    count = 0
    children_idx = children(g, k)
    for i in children_idx
        is_num(b, i) && continue                     # Don't contract a number valued argument
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
        _, w, _ = IntervalContractors.mul_rev(interval(b, k), interval(b, i), Intv(tmul))
        isempty(w) && (return false)
        b[i] = MC{N,T}(w)
    end
    return true
end

for (f, fc, F) in ((^, POW, IntervalContractors.power_rev),
                   (/, DIV, IntervalContractors.div_rev))
    @eval function rprop!(t::Relax, v::Val{$fc}, g::DAT, b::RelaxCache{V,N,T}, k) where {V,N,T<:RelaxTag}
        is_num(b, k) && (return true)
        x = child(g, 1, k)
        y = child(g, 2, k)

        if xset_ynum(b, x, y)
            z, u, v = ($F)(interval(b, k), interval(b, x), num(b, y))
        elseif xnum_yset(b, x, y)
            z, u, v = ($F)(interval(b, k), num(b, x), interval(b, y))
        else
            z, u, v = ($F)(interval(b, k), interval(b, x), interval(b, y))
        end
        if !is_num(b, x)
            isempty(u) && (return false)
            b[x] = MC{N,T}(u)
        end
        if !is_num(b, y)
            isempty(v) && (return false)
            b[y] = MC{N,T}(v)
        end
        return true
    end
end

rprop!(t::Relax, v::Val{USER}, g::DAT, b::RelaxCache, k::Int) = true
rprop!(t::Relax, v::Val{USERN}, g::DAT, b::RelaxCache, k::Int) = true

for (fc, F) in ((MINUS, IntervalContractors.minus_rev),
                (SQRT, IntervalContractors.sqrt_rev),
                (ABS, IntervalContractors.abs_rev),
                (EXP, IntervalContractors.exp_rev),
                (EXP2, IntervalContractors.exp2_rev),
                (EXP10, IntervalContractors.exp10_rev),
                (EXPM1, IntervalContractors.expm1_rev),
                (LOG, IntervalContractors.log_rev),
                (LOG2, IntervalContractors.log2_rev),
                (LOG10, IntervalContractors.log10_rev),
                (LOG1P, IntervalContractors.log1p_rev),
                (SIN, IntervalContractors.sin_rev),
                (COS, IntervalContractors.cos_rev),
                (TAN, IntervalContractors.tan_rev),
                (ASIN, IntervalContractors.asin_rev),
                (ACOS, IntervalContractors.acos_rev),
                (ATAN, IntervalContractors.atan_rev),
                (SINH, IntervalContractors.sinh_rev),
                (COSH, IntervalContractors.cosh_rev),
                (TANH, IntervalContractors.tanh_rev),
                (ASINH, IntervalContractors.asinh_rev),
                (ACOSH, IntervalContractors.acosh_rev),
                (ATANH, IntervalContractors.atanh_rev),
                )
    @eval function rprop!(t::Relax, v::Val{$fc}, g::DAT, b::RelaxCache{V,N,T}, k) where {V,N,T<:RelaxTag}
        is_num(b, k) && (return true)
        x = child(g, 1, k)
        z, u = ($F)(interval(b, k), interval(b, x))
        if !is_num(b, x)
            isempty(u) && (return false)
            b[x] = MC{N,T}(u)
        end
        return true
    end
end
