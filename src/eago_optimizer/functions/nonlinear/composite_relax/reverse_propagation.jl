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

function r_init!(t::Relax, tree::DAT, cache::RelaxCache{V,N,T}) where {V,N,T<:RelaxTag}
    if !is_num(cache, 1)
        cache[1] = set(cache, 1) ∩ tree.sink_bnd
    end
    return !isempty(set(cache, 1))
end

function rprop!(t::RelaxCacheAttribute, v::Variable, tree::DAT, cache::RelaxCache{V,N,T}, node_index::Int) where {V,N,T<:RelaxTag}
    variable_index = first_index(tree, node_index)
    variable_value = val(cache, variable_index)
    variable_lower_bound = lbd(cache, variable_index)
    variable_upper_bound = ubd(cache, variable_index)
    node_interval = interval(cache, node_index)
    if variable_lower_bound == variable_upper_bound
        cache[node_index] = variable_value
    else
        variable_node = varset(MC{N,T}, rev_sparsity(tree, variable_index, node_index), variable_value, variable_value, variable_lower_bound, variable_upper_bound)
        variable_node = variable_node ∩ node_interval
        cache[node_index] = variable_node
        cache.ic.v.lower_variable_bounds[variable_index] = variable_node.Intv.lo
        cache.ic.v.upper_variable_bounds[variable_index] = variable_node.Intv.hi
    end
    return !isempty(set(cache, node_index))
end

function rprop!(t::Relax, ex::Subexpression, tree::DAT, cache::RelaxCache{V,N,T}, node_index::Int) where {V,N,T<:RelaxTag}
    first_child_index = first_index(tree, node_index)
    subexpression_index = dependent_subexpression_index(tree, first_child_index)
    if subexpression_is_num(cache, subexpression_index)
        cache[node_index] = subexpression_num(cache, subexpression_index)
    else
        cache[node_index] = subexpression_set(cache, subexpression_index)
    end
    return true
end

for F in (PLUS, MULT)
    @eval function rprop!(t::Relax, v::Val{$F}, tree::DAT, cache::RelaxCache{V,N,T}, node_index::Int) where {V,N,T<:RelaxTag}
        is_binary(tree, node_index) ? rprop_2!(Relax(), Val($F), tree, cache, node_index) : rprop_n!(Relax(), Val($F), tree, cache, node_index)
    end
end

const MAX_ASSOCIATIVE_REVERSE = 6

"""
$(FUNCTIONNAME)

Updates storage tapes with reverse evaluation of node representing `n = x + y` which updates x and y.
"""
function rprop_2!(t::Relax, v::Val{PLUS}, tree::DAT, cache::RelaxCache{V,N,T}, node_index::Int) where {V,N,T<:RelaxTag}
    is_num(cache, node_index) && return true
    first_child_index = child(tree, 1, node_index)
    second_child_index = child(tree, 2, node_index)
    node_interval = interval(cache, node_index)
    old_first_child_interval = interval(cache, first_child_index)
    old_second_child_interval = interval(cache, second_child_index)

    if xset_ynum(cache, first_child_index, second_child_index)
        _, new_first_child_interval, new_second_child_interval = IntervalContractors.plus_rev(node_interval, old_first_child_interval, num(cache, second_child_index))
    elseif xnum_yset(cache, first_child_index, second_child_index)
        _, new_first_child_interval, new_second_child_interval = IntervalContractors.plus_rev(node_interval, num(cache, first_child_index), old_second_child_interval)
    else
        _, new_first_child_interval, new_second_child_interval = IntervalContractors.plus_rev(node_interval, old_first_child_interval, old_second_child_interval)
    end
        if !is_num(cache, first_child_index)
            isempty(new_first_child_interval) && return false
            cache[first_child_index] = MC{N,T}(new_first_child_interval)
        end
        if !is_num(cache, second_child_index)
            isempty(new_second_child_interval) && return false
            cache[second_child_index] = MC{N,T}(new_second_child_interval)
        end
    return true
end

"""
$(FUNCTIONNAME)

Updates storage tapes with reverse evaluation of node representing `n = +(x,y,z...)` which updates x, y, z and so on.
"""
function rprop_n!(t::Relax, v::Val{PLUS}, tree::DAT, cache::RelaxCache{V,N,T}, node_index::Int) where {V,N,T<:RelaxTag}
    # Outer loop makes a temporary sum (minus one argument)
    # A reverse is then computed with respect to this argument
    count = 0
    children_indices = children(tree, node_index)
    node_interval = interval(cache, node_index)
    for current_child_index in children_indices
        is_num(cache, current_child_index) && continue                     # Don't contract a number valued argument
        (count >= MAX_ASSOCIATIVE_REVERSE) && break
        old_current_child_interval = interval(cache, current_child_index)
        temp_sum = zero(MC{N,T})
        count += 1
        for other_child_index in children_indices
            if other_child_index != current_child_index
                if is_num(cache, other_child_index)
                    temp_sum += num(cache, other_child_index)
                else
                    temp_sum += set(cache, other_child_index)
                end
            end
        end
        _, new_current_child_interval, _ = IntervalContractors.plus_rev(node_interval, old_current_child_interval, Intv(temp_sum))
        isempty(new_current_child_interval) && (return false)
        cache[current_child_index] = MC{N,T}(new_current_child_interval)
    end
    return true
end

"""
$(FUNCTIONNAME)

Updates storage tapes with reverse evaluation of node representing `n = x * y` which updates x and y.
"""
function rprop_2!(t::Relax, v::Val{MULT}, tree::DAT, cache::RelaxCache{V,N,T}, node_index::Int) where {V,N,T<:RelaxTag}

    is_num(cache, node_index) && (return true)
    first_child_index = child(tree, 1, node_index)
    second_child_index = child(tree, 2, node_index)
    node_interval = interval(cache, node_index)
    old_first_child_interval = interval(cache, first_child_index)
    old_second_child_interval = interval(cache, second_child_index)

    if xset_ynum(cache, first_child_index, second_child_index)
        _, new_first_child_interval, new_second_child_interval = IntervalContractors.mul_rev(node_interval, old_first_child_interval, num(cache, second_child_index))
    elseif xnum_yset(cache, first_child_index, second_child_index)
        _, new_first_child_interval, new_second_child_interval = IntervalContractors.mul_rev(node_interval, num(cache, first_child_index), old_second_child_interval)
    else
        _, new_first_child_interval, new_second_child_interval = IntervalContractors.mul_rev(node_interval, old_first_child_interval, old_second_child_interval)
    end

    if !is_num(cache, first_child_index)
        isempty(new_first_child_interval) && (return false)
        cache[first_child_index] = MC{N,T}(new_first_child_interval)
    end
    if !is_num(cache, second_child_index)
        isempty(new_second_child_interval) && (return false)
        cache[second_child_index] = MC{N,T}(new_second_child_interval)
    end
    return true
end

"""
$(FUNCTIONNAME)

Updates storage tapes with reverse evaluation of node representing `n = *(x,y,z...)` which updates x, y, z and so on.
"""
function rprop_n!(t::Relax, v::Val{MULT}, tree::DAT, cache::RelaxCache{V,N,T}, node_index::Int) where {V,N,T<:RelaxTag}
    # Outer loop makes a temporary product (minus one argument)
    # A reverse is then computed with respect to this argument
    count = 0
    children_indices = children(tree, node_index)
    node_interval = interval(cache, node_index)
    for current_child_index in children_indices
        is_num(cache, current_child_index) && continue                     # Don't contract a number valued argument
        (count >= MAX_ASSOCIATIVE_REVERSE) && break
        old_current_child_interval = interval(cache, current_child_index)
        temp_prod = one(MC{N,T})
        count += 1
        for other_child_index in children_indices
            if current_child_index != other_child_index
                if is_num(cache, other_child_index)
                    temp_prod *= num(cache, other_child_index)
                else
                    temp_prod *= set(cache, other_child_index)
                end
            end
        end
        _, new_current_child_interval, _ = IntervalContractors.mul_rev(node_interval, old_current_child_interval, Intv(temp_prod))
        isempty(new_current_child_interval) && (return false)
        cache[current_child_index] = MC{N,T}(new_current_child_interval)
    end
    return true
end

function rprop!(t::Relax, v::Val{MINUS}, tree::DAT, cache::RelaxCache{V,N,T}, node_index::Int) where {V,N,T<:RelaxTag}
    is_num(cache, node_index) && (return true)
    if is_unary(tree, node_index)
        first_child_index = child(tree, 1, node_index)
        node_interval = interval(cache, node_index)
        old_first_child_interval = interval(cache, first_child_index)
        _, new_first_child_interval = IntervalContractors.minus_rev(node_interval, old_first_child_interval)
        if !is_num(cache, first_child_index)
            isempty(new_first_child_interval) && (return false)
            cache[first_child_index] = MC{N,T}(new_first_child_interval)
        end
    else
        first_child_index = child(tree, 1, node_index)
        second_child_index = child(tree, 2, node_index)
        node_interval = interval(cache, node_index)
        old_first_child_interval = interval(cache, first_child_index)
        old_second_child_interval = interval(cache, second_child_index)

        if xset_ynum(cache, first_child_index, second_child_index)
            _, new_first_child_interval, new_second_child_interval = IntervalContractors.minus_rev(node_interval, old_first_child_interval, num(cache, second_child_index))
        elseif xnum_yset(cache, first_child_index, second_child_index)
            _, new_first_child_interval, new_second_child_interval = IntervalContractors.minus_rev(node_interval, num(cache, first_child_index), old_second_child_interval)
        else
            _, new_first_child_interval, new_second_child_interval = IntervalContractors.minus_rev(node_interval, old_first_child_interval, old_second_child_interval)
        end
        if !is_num(cache, first_child_index)
            isempty(new_first_child_interval) && (return false)
            cache[first_child_index] = MC{N,T}(new_first_child_interval)
        end
        if !is_num(cache, second_child_index)
            isempty(new_second_child_interval) && (return false)
            cache[second_child_index] = MC{N,T}(new_second_child_interval)
        end
    end
    return true
end

for (f, fc, F) in ((^, POW, IntervalContractors.power_rev),
                   (/, DIV, IntervalContractors.div_rev))
    @eval function rprop!(t::Relax, v::Val{$fc}, tree::DAT, cache::RelaxCache{V,N,T}, node_index::Int) where {V,N,T<:RelaxTag}
        is_num(cache, node_index) && (return true)
        first_child_index = child(tree, 1, node_index)
        second_child_index = child(tree, 2, node_index)
        node_interval = interval(cache, node_index)
        old_first_child_interval = interval(cache, first_child_index)
        old_second_child_interval = interval(cache, second_child_index)

        if xset_ynum(cache, first_child_index, second_child_index)
            _, new_first_child_interval, new_second_child_interval = ($F)(node_interval, old_first_child_interval, num(cache, second_child_index))
        elseif xnum_yset(cache, first_child_index, second_child_index)
            _, new_first_child_interval, new_second_child_interval = ($F)(node_interval, num(cache, first_child_index), old_second_child_interval)
        else
            _, new_first_child_interval, new_second_child_interval = ($F)(node_interval, old_first_child_interval, old_second_child_interval)
        end
        if !is_num(cache, first_child_index)
            isempty(new_first_child_interval) && (return false)
            cache[first_child_index] = MC{N,T}(new_first_child_interval)
        end
        if !is_num(cache, second_child_index)
            isempty(new_second_child_interval) && (return false)
            cache[second_child_index] = MC{N,T}(new_second_child_interval)
        end
        return true
    end
end

rprop!(t::Relax, v::Val{USER}, tree::DAT, cache::RelaxCache, k::Int) = true
rprop!(t::Relax, v::Val{USERN}, tree::DAT, cache::RelaxCache, k::Int) = true

for (fc, F) in ((SQRT, IntervalContractors.sqrt_rev),
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
    @eval function rprop!(t::Relax, v::Val{$fc}, tree::DAT, cache::RelaxCache{V,N,T}, node_index::Int) where {V,N,T<:RelaxTag}
        is_num(cache, node_index) && (return true)
        first_child_index = child(tree, 1, node_index)
        node_interval = interval(cache, node_index)
        old_first_child_interval = interval(cache, first_child_index)
        _, new_first_child_interval = ($F)(node_interval, old_first_child_interval)
        if !is_num(cache, first_child_index)
            isempty(new_first_child_interval) && (return false)
            cache[first_child_index] = MC{N,T}(new_first_child_interval)
        end
        return true
    end
end
