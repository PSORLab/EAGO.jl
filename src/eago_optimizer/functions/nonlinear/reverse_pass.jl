# Copyright (c) 2018: Matthew Wilhelm & Matthew Stuber.
# This work is licensed under the Creative Commons Attribution-NonCommercial-
# ShareAlike 4.0 International License. To view a copy of this license, visit
# http://creativecommons.org/licenses/by-nc-sa/4.0/ or send a letter to Creative
# Commons, PO Box 1866, Mountain View, CA 94042, USA.
#############################################################################
# EAGO
# A development environment for robust and global optimization
# See https://github.com/PSORLab/EAGO.jl
#############################################################################
# src/eago_optimizer/functions/nonlinear/reverse_pass.jl
# Functions used to compute reverse pass of nonlinear functions.
#############################################################################

# maximum number to perform reverse operation on associative term by summing
# and evaluating pairs remaining terms not reversed
const MAX_ASSOCIATIVE_REVERSE = 6
const REVERSE_DEBUG = false

"""
$(FUNCTIONNAME)

Updates storage tapes with reverse evalution of node representing `n = x + y` which updates x & y.
"""
function reverse_plus_binary!(k::Int64, children_arr::Vector{Int64}, children_idx::UnitRange{Int64},
                              numvalued::Vector{Bool}, numberstorage::Vector{Float64}, setstorage::Vector{MC{N,T}},
                              x::Vector{Float64}, lbd::Vector{Float64}, ubd::Vector{Float64}, sparsity::Vector{Int},
                              subgrad_tol::Float64, is_post::Bool) where {N, T<:RelaxTag}

    REVERSE_DEBUG && println("--- start reverse plus binary ---")

    # extract values for k
    argk_is_number =  numvalued[k]
    if !argk_is_number
        setk =  setstorage[k]
    end

    # get row indices
    idx1 = first(children_idx)
    idx2 = last(children_idx)

    # extract values for argument 1
    arg1_index = children_arr[idx1]
    arg1_is_number = numvalued[arg1_index]
    if arg1_is_number
        set1 = zero(MC{N,T})
        num1 = numberstorage[arg1_index]
    else
        num1 = 0.0
        set1 = setstorage[arg1_index]
    end

    # extract values for argument 2
    arg2_index =  children_arr[idx2]
    arg2_is_number =  numvalued[arg2_index]
    if arg2_is_number
        num2 = numberstorage[arg2_index]
        set2 = zero(MC{N,T})
    else
        set2 = setstorage[arg2_index]
        num2 = 0.0
    end

    if !arg1_is_number && arg2_is_number
        c, a, b = plus_rev(setk, set1, num2)

    elseif arg1_is_number && !arg2_is_number
        c, a, b = plus_rev(setk, num1, set2)

    else
        c, a, b = plus_rev(setk, set1, set2)
    end

    if REVERSE_DEBUG
        println("val out = $(c)")
        println("arg1 out = $(a)")
        println("arg2 out = $(b)")
    end

    if !arg1_is_number
        if isempty(a)
            return false
        elseif isnan(a)
            setstorage[arg1_index] = MC{N,T}(a.Intv)
        else
            setstorage[arg1_index] = is_post ? set_value_post(x, a, lbd, ubd, sparsity, subgrad_tol) : a
        end
        REVERSE_DEBUG && println("setstorage[arg1_index] = $(setstorage[arg1_index])")
    end

    if !arg2_is_number
        if isempty(b)
            return false
        elseif isnan(b)
            setstorage[arg2_index] = MC{N,T}(b.Intv)
        else
            setstorage[arg2_index] = is_post ? set_value_post(x, b, lbd, ubd, sparsity, subgrad_tol) : b
        end
        REVERSE_DEBUG && println("setstorage[arg2_index] = $(setstorage[arg2_index])")
    end

    return true
end

"""
$(FUNCTIONNAME)

Updates storage tapes with reverse evalution of node representing `n = +(x,y,z...)` which updates x, y, z and so on.
"""
function reverse_plus_narity!(k::Int64, children_arr::Vector{Int64}, children_idx::UnitRange{Int64},
                              numvalued::Vector{Bool}, numberstorage::Vector{Float64}, setstorage::Vector{MC{N,T}},
                              x::Vector{Float64}, lbd::Vector{Float64}, ubd::Vector{Float64}, sparsity::Vector{Int},
                              subgrad_tol::Float64, is_post::Bool) where {N, T<:RelaxTag}

    REVERSE_DEBUG && println("--- start reverse plus narity ---")

    continue_flag = true
    argk_is_number = numvalued[k]
    if !argk_is_number
        setk =  setstorage[k]
    end

    # out loops makes a temporary sum (minus one argument)
    # a reverse is then compute with respect to this argument
    active_count_number = 0
    for idx in children_idx
        active_idx = children_arr[idx]

        # don't contract a number valued argument
        active_arg_is_number = numvalued[active_idx]
        active_arg_is_number && continue

        if active_count_number >= MAX_ASSOCIATIVE_REVERSE
            break
        end

        tmp_sum = zero(MC{N,T})
        active_count_number += 1
        for nidx in children_idx
            inactive_idx =  children_arr[nidx]
            if inactive_idx != active_idx
                if  numvalued[inactive_idx]
                    tmp_sum += numberstorage[inactive_idx]
                else
                    tmp_sum += setstorage[inactive_idx]
                end
            end
        end

        active_set = setstorage[active_idx]
        c, a, b = plus_rev(setk, active_set, tmp_sum)

        if REVERSE_DEBUG
            println("val out = $(c)")
            println("arg1 out = $(a)")
            println("arg2 out = $(b)")
        end

        if isempty(a)
            return false
        elseif isnan(a)
            setstorage[active_idx] = MC{N,T}(a.Intv)
        else
            setstorage[active_idx] = is_post ? set_value_post(x, a, lbd, ubd, sparsity, subgrad_tol) : a
        end
    end

    return true
end

"""
$(FUNCTIONNAME)

Updates storage tapes with reverse evalution of node representing `n = x * y` which updates x & y.
"""
function reverse_multiply_binary!(k::Int64, children_arr::Vector{Int64}, children_idx::UnitRange{Int64},
                                  numvalued::Vector{Bool}, numberstorage::Vector{Float64}, setstorage::Vector{MC{N,T}},
                                  x::Vector{Float64}, lbd::Vector{Float64}, ubd::Vector{Float64}, sparsity::Vector{Int},
                                  subgrad_tol::Float64, is_post::Bool) where {N, T<:RelaxTag}

    REVERSE_DEBUG && println(" is_post = $is_post")

    REVERSE_DEBUG && println("--- start reverse mult binary ---")

    # extract values for k
    argk_is_number = numvalued[k]
    if !argk_is_number
        setk = setstorage[k]
        REVERSE_DEBUG && println("setk = $setk")
    end

    # get row indices
    idx1 = first(children_idx)
    idx2 = last(children_idx)

    # extract values for argument 1
    arg1_index = children_arr[idx1]
    arg1_is_number = numvalued[arg1_index]
    if arg1_is_number
        set1 = zero(MC{N,T})
        num1 =  numberstorage[arg1_index]
    else
        num1 = 0.0
        set1 = setstorage[arg1_index]
    end

    # extract values for argument 2
    arg2_index =  children_arr[idx2]
    arg2_is_number = numvalued[arg2_index]
    if arg2_is_number
        num2 = numberstorage[arg2_index]
        set2 = zero(MC{N,T})
    else
        set2 = setstorage[arg2_index]
        num2 = 0.0
    end

    if !arg1_is_number && arg2_is_number
        c, a, b = mult_rev(setk, set1, num2)

    elseif arg1_is_number && !arg2_is_number
        c, a, b = mult_rev(setk, num1, set2)

    else
        c, a, b = mult_rev(setk, set1, set2)
    end

    if REVERSE_DEBUG
        println("val out = $(c)")
        println("arg1 out = $(a)")
        println("arg2 out = $(b)")
    end

    if !arg1_is_number
        if isempty(a)
            return false
        elseif isnan(a)
            setstorage[arg1_index] = MC{N,T}(a.Intv)
        else
            setstorage[arg1_index] = is_post ? set_value_post(x, a, lbd, ubd, sparsity, subgrad_tol) : a
        end
        REVERSE_DEBUG && println("setstorage[arg1_index] = $(setstorage[arg1_index])")
    end

    if !arg2_is_number
        if isempty(b)
            return false
        elseif isnan(b)
            setstorage[arg2_index] = MC{N,T}(b.Intv)
        else
            setstorage[arg2_index] = is_post ? set_value_post(x, b, lbd, ubd, sparsity, subgrad_tol) : b
        end
        REVERSE_DEBUG && println("setstorage[arg2_index] = $(setstorage[arg2_index])")
    end

    return true
end

"""
$(FUNCTIONNAME)

Updates storage tapes with reverse evalution of node representing `n = *(x,y,z...)` which updates x, y, z and so on.
"""
function reverse_multiply_narity!(k::Int64, children_arr::Vector{Int64}, children_idx::UnitRange{Int64},
                              numvalued::Vector{Bool}, numberstorage::Vector{Float64}, setstorage::Vector{MC{N,T}},
                              x::Vector{Float64}, lbd::Vector{Float64}, ubd::Vector{Float64}, sparsity::Vector{Int},
                              subgrad_tol::Float64, is_post::Bool) where {N, T<:RelaxTag}

    REVERSE_DEBUG && println("--- start reverse mult narity ---")
    continue_flag = true
    argk_is_number = numvalued[k]
    if !argk_is_number
        setk = setstorage[k]
    end

    # out loops makes a temporary sum (minus one argument)
    # a reverse is then compute with respect to this argument
    active_count_number = 0
    for idx in children_idx
        active_idx = children_arr[idx]

        # don't contract a number valued argument
        active_arg_is_number = numvalued[active_idx]
        active_arg_is_number && continue

        if active_count_number >= MAX_ASSOCIATIVE_REVERSE
            break
        end

        tmp_mul = one(MC{N,T})
        active_count_number += 1
        for nidx in children_idx
            inactive_idx = children_arr[nidx]
            if inactive_idx != active_idx
                if  numvalued[inactive_idx]
                    tmp_mul *=  numberstorage[inactive_idx]
                else
                    tmp_mul *=  setstorage[inactive_idx]
                end
            end
        end

        active_set = setstorage[active_idx]
        c, a, b = mult_rev(setk, active_set, tmp_mul)

        if REVERSE_DEBUG
            println("val out = $(c)")
            println("arg1 out = $(a)")
            println("arg2 out = $(b)")
        end

        if isempty(a)
            return false
        elseif isnan(a)
            setstorage[active_idx] = MC{N,T}(a.Intv)
        else
            setstorage[active_idx] = is_post ? set_value_post(x, a, lbd, ubd, sparsity, subgrad_tol) : a
        end
    end

    return true
end

function reverse_minus!(k::Int64, children_arr::Vector{Int64}, children_idx::UnitRange{Int64},
                        numvalued::Vector{Bool}, numberstorage::Vector{Float64}, setstorage::Vector{MC{N,T}},
                        x::Vector{Float64}, lbd::Vector{Float64}, ubd::Vector{Float64}, sparsity::Vector{Int},
                        subgrad_tol::Float64, is_post::Bool) where {N, T<:RelaxTag}

    REVERSE_DEBUG && println("--- start reverse minus ---")

    #argk_index =  children_arr[k]
    argk_is_number = numvalued[k]
    if !argk_is_number
        setk = setstorage[k]
    end

    # don't perform a reverse pass if the output was a number
    if argk_is_number
        return true
    end

    # get row indices
    idx1 = first(children_idx)
    idx2 = last(children_idx)

    # extract values for argument 1
    arg1_index = children_arr[idx1]
    arg1_is_number = numvalued[arg1_index]
    if arg1_is_number
        set1 = zero(MC{N,T})
        num1 = numberstorage[arg1_index]
    else
        num1 = 0.0
        set1 = setstorage[arg1_index]
    end

    # extract values for argument 2
    arg2_index = children_arr[idx2]
    arg2_is_number = numvalued[arg2_index]
    if arg2_is_number
        num2 = numberstorage[arg2_index]
        set2 = zero(MC{N,T})
    else
        set2 = setstorage[arg2_index]
        num2 = 0.0
    end

    if !arg1_is_number && arg2_is_number
        c, a, b = minus_rev(setk, set1, num2)

    elseif arg1_is_number && !arg2_is_number
        c, a, b = minus_rev(setk, num1, set2)

    else
        c, a, b = minus_rev(setk, set1, set2)
    end

    if REVERSE_DEBUG
        println("val out = $(c)")
        println("arg1 out = $(a)")
        println("arg2 out = $(b)")
    end

    if !arg1_is_number
        if isempty(a)
            return false
        elseif isnan(a)
            a = MC{N,T}(a.Intv)
            setstorage[arg1_index] = a
        else
            setstorage[arg1_index] = is_post ? set_value_post(x, a, lbd, ubd, sparsity, subgrad_tol) : a
        end
        REVERSE_DEBUG && println("setstorage[arg1_index] = $(setstorage[arg1_index])")
    end

    if !arg2_is_number
        if isempty(b)
            return false
        elseif isnan(b)
            b = MC{N,T}(b.Intv)
            setstorage[arg2_index] = b
        else
            setstorage[arg2_index] = is_post ? set_value_post(x, b, lbd, ubd, sparsity, subgrad_tol) : b
        end
        REVERSE_DEBUG && println("setstorage[arg2_index] = $(setstorage[arg2_index])")
    end

    return true
end

function reverse_power!(k::Int64, children_arr::Vector{Int64}, children_idx::UnitRange{Int64},
                        numvalued::Vector{Bool}, numberstorage::Vector{Float64}, setstorage::Vector{MC{N,T}},
                        x::Vector{Float64}, lbd::Vector{Float64}, ubd::Vector{Float64}, sparsity::Vector{Int},
                        subgrad_tol::Float64, is_post::Bool) where {N, T<:RelaxTag}

    REVERSE_DEBUG && println("--- start reverse power ---")

    # extract values for k
    argk_is_number =  numvalued[k]
    if !argk_is_number
        setk = setstorage[k]
    end

    # don't perform a reverse pass if the output was a number
    if argk_is_number
        return true
    end

    # get row indices
    idx1 = first(children_idx)
    idx2 = last(children_idx)

    # extract values for argument 1
    arg1_index = children_arr[idx1]
    arg1_is_number = numvalued[arg1_index]
    if arg1_is_number
        set1 = zero(MC{N,T})
        num1 = numberstorage[arg1_index]
    else
        num1 = 0.0
        set1 = setstorage[arg1_index]
    end

    # extract values for argument 2
    arg2_index = children_arr[idx2]
    arg2_is_number = numvalued[arg2_index]
    if arg2_is_number
        num2 = numberstorage[arg2_index]
        set2 = zero(MC{N,T})
    else
        set2 = setstorage[arg2_index]
        num2 = 0.0
    end

    if !arg1_is_number && arg2_is_number
        c, a, b = power_rev(setk, set1, num2)

    elseif arg1_is_number && !arg2_is_number
        c, a, b = power_rev(setk, num1, set2)

    else
        c, a, b = power_rev(setk, set1, set2)
    end

    if REVERSE_DEBUG
        println("val out = $(c)")
        println("arg1 out = $(a)")
        println("arg2 out = $(b)")
    end

    if !arg1_is_number
        if isempty(a)
            return false
        elseif isnan(a)
            setstorage[arg1_index] = MC{N,T}(a.Intv)
        else
            setstorage[arg1_index] = is_post ? set_value_post(x, a, lbd, ubd) : a
        end
    end

    if !arg2_is_number
        if isempty(b)
            return false
        elseif isnan(b)
            setstorage[arg2_index] = MC{N,T}(b.Intv)
        else
            setstorage[arg2_index] = is_post ? set_value_post(x, b, lbd, ubd, sparsity, subgrad_tol) : b
        end
    end

    return true
end

function reverse_divide!(k::Int64, children_arr::Vector{Int64}, children_idx::UnitRange{Int64},
                        numvalued::Vector{Bool}, numberstorage::Vector{Float64}, setstorage::Vector{MC{N,T}},
                        x::Vector{Float64}, lbd::Vector{Float64}, ubd::Vector{Float64}, sparsity::Vector{Int},
                        subgrad_tol::Float64, is_post::Bool) where {N, T<:RelaxTag}

    REVERSE_DEBUG && println("--- start reverse divide ---")

    # extract values for k
    argk_is_number =  numvalued[k]
    if !argk_is_number
        setk = setstorage[k]
    end

    # don't perform a reverse pass if the output was a number
    if argk_is_number
        return true
    end

    # get row indices
    idx1 = first(children_idx)
    idx2 = last(children_idx)

    # extract values for argument 1
    arg1_index = children_arr[idx1]
    arg1_is_number = numvalued[arg1_index]
    if arg1_is_number
        set1 = zero(MC{N,T})
        num1 = numberstorage[arg1_index]
    else
        num1 = 0.0
        set1 = setstorage[arg1_index]
    end

    # extract values for argument 2
    arg2_index = children_arr[idx2]
    arg2_is_number = numvalued[arg2_index]
    if arg2_is_number
        num2 = numberstorage[arg2_index]
        set2 = zero(MC{N,T})
    else
        set2 =  setstorage[arg2_index]
        num2 = 0.0
    end

    if !arg1_is_number && arg2_is_number
        c, a, b = div_rev(setk, set1, num2)

    elseif arg1_is_number && !arg2_is_number
        c, a, b = div_rev(setk, num1, set2)

    else
        c, a, b = div_rev(setk, set1, set2)
    end

    if REVERSE_DEBUG
        println("val out = $(c)")
        println("arg1 out = $(a)")
        println("arg2 out = $(b)")
    end

    if !arg1_is_number
        if isempty(a)
            return false
        elseif isnan(a)
            setstorage[arg1_index] = MC{N,T}(a.Intv)
        else
            setstorage[arg1_index] = is_post ? set_value_post(x, a, lbd, ubd, sparsity, subgrad_tol) : a
        end
    end

    if !arg2_is_number
        if isempty(b)
            return false
        elseif isnan(b)
            setstorage[arg2_index] = MC{N,T}(b.Intv)
        else
            setstorage[arg2_index] = is_post ? set_value_post(x, b, lbd, ubd, sparsity, subgrad_tol) : b
        end
    end

    return true
end

function reverse_univariate!(k::Int64, op::Int64, arg_indx::Int64, setstorage::Vector{MC{N,T}}, x::Vector{Float64},
                             lbd::Vector{Float64}, ubd::Vector{Float64}, sparsity::Vector{Int},
                             subgrad_tol::Float64, is_post::Bool) where {N, T<:RelaxTag}

    REVERSE_DEBUG && println("--- start reverse reverse_univariate ---")
    valset = setstorage[k]
    argset = setstorage[arg_indx]

    REVERSE_DEBUG && println("val set = $(valset)")
    REVERSE_DEBUG && println("arg set = $(argset)")

    a, b = eval_univariate_set_reverse(op, valset, argset)

    REVERSE_DEBUG && println("val out = $(a)")
    REVERSE_DEBUG && println("arg out = $(b)")

    if isempty(b)
        return false
    elseif isnan(b)
        setstorage[arg_indx] = MC{N,T}(b.Intv)
    else
        setstorage[arg_indx] = is_post ? set_value_post(x, b, lbd, ubd, sparsity, subgrad_tol) : b
    end

    return true
end

function reverse_set_subexpression!(k::Int64, op::Int64, subexpressions::Vector{NonlinearExpression},
                                    numvalued::Vector{Bool}, numberstorage::Vector{Float64},
                                    setstorage::Vector{MC{N,T}}, cv_buffer::Vector{Float64},
                                    cc_buffer::Vector{Float64}, func_sparsity::Vector{Int64}) where {N, T<:RelaxTag}

    is_number = subexpression_isnum[nod.index]
    if !is_number
         subexpr_values_set[nod.index] = setstorage[k]
    end

    subexpression = subexpressions[op]

    isa_number = subexpression.is_number[1]
    if !isa_number
        copy_subexpression_value!(k, op, setstorage, subexpression, cv_grad_buffer, cc_grad_buffer)
    end
    numvalued[k] = isa_number

    return nothing
end

"""
$(TYPEDSIGNATURES)

Performs a reverse McCormick/interval pass. If a NaN value is computed for the McCormick relaxation then the
routine defaults to the interval value instead.

There is a tacit assumption here that an abstract tree structure is used for propagation. If a cse structure is
eventually used then each reverse_xxx function will need to keep track of each variables state.
"""
function reverse_pass_kernel!(nd::Vector{JuMP.NodeData}, adj::SparseMatrixCSC{Bool,Int64}, x::Vector{Float64},
                              lbd::Vector{Float64}, ubd::Vector{Float64}, sparsity::Vector{Int},
                              setstorage::Vector{MC{N,T}}, subgrad_tol::Float64,
                              numberstorage::Vector{Float64}, numvalued::Vector{Bool},
                              is_post_input::Bool) where {N, T<:RelaxTag}

    children_arr = rowvals(adj)
    continue_flag = true
    is_post = is_post_input

    for k = 1:length(nd)

        REVERSE_DEBUG && println(" ")
        nod = nd[k]
        ntype = nod.nodetype
        nvalued = numvalued[k]

        if ntype == JuMP._Derivatives.VALUE      || ntype == JuMP._Derivatives.LOGIC     ||
           ntype == JuMP._Derivatives.COMPARISON || ntype == JuMP._Derivatives.PARAMETER ||
           ntype == JuMP._Derivatives.EXTRA      || ntype == JuMP._Derivatives.SUBEXPRESSION
           continue

        elseif nod.nodetype == JuMP._Derivatives.VARIABLE
            op = nod.index
            REVERSE_DEBUG && println("--- start reverse reverse variable ---")
            variable_interval = setstorage[k].Intv
            lower_interval = variable_interval.lo
            upper_interval = variable_interval.hi

            # update
            prior_x = x[op]
            if lower_interval > prior_x
                x[op] = 0.5*(lower_interval + upper_interval)
                is_post = false
            end
            if upper_interval < prior_x
                x[op] = 0.5*(lower_interval + upper_interval)
                is_post = false
            end
            lbd[op] = lower_interval
            ubd[op] = upper_interval

            REVERSE_DEBUG && println("variable_rev[$op][$k] at k = $k -> $(setstorage[k])")

        elseif nvalued
            continue

        elseif nod.nodetype == JuMP._Derivatives.CALL
            op = nod.index
            parent_index = nod.parent
            children_idx = nzrange(adj, k)
            parent_value = setstorage[k]
            n_children = length(children_idx)

            # SKIPS USER DEFINE OPERATORS NOT BRIDGED INTO JuMP Tree Representation
            if op >= JuMP._Derivatives.USER_OPERATOR_ID_START
                continue

            # :+
            elseif op === 1
                REVERSE_DEBUG && println("plus_rev[$n_children][$k]   at k = $k -> $(setstorage[k])")
                if n_children === 2
                    continue_flag &= reverse_plus_binary!(k, children_arr, children_idx, numvalued, numberstorage,
                                                          setstorage, x, lbd, ubd, sparsity, subgrad_tol, is_post)
                else
                    continue_flag &= reverse_plus_narity!(k, children_arr, children_idx, numvalued, numberstorage,
                                                          setstorage, x, lbd, ubd, sparsity, subgrad_tol, is_post)
                end
                REVERSE_DEBUG && !continue_flag && println("Infeasible node encountered.")

            # :-
            elseif op === 2
                REVERSE_DEBUG && println("minus_rev[$k]        at k = $k -> $(setstorage[k])")
                continue_flag &= reverse_minus!(k, children_arr, children_idx, numvalued, numberstorage,
                                                setstorage, x, lbd, ubd, sparsity, subgrad_tol, is_post)
                REVERSE_DEBUG && !continue_flag && println("Infeasible node encountered.")

            elseif op === 3 # :*
                REVERSE_DEBUG && println("mult_rev[$n_children][$k]     at k = $k -> $(setstorage[k])")
                if n_children === 2
                    continue_flag &= reverse_multiply_binary!(k, children_arr, children_idx, numvalued, numberstorage,
                                                              setstorage, x, lbd, ubd, sparsity, subgrad_tol, is_post)
                else
                    continue_flag &= reverse_multiply_narity!(k, children_arr, children_idx, numvalued, numberstorage,
                                                              setstorage, x, lbd, ubd, sparsity, subgrad_tol, is_post)
                end
                REVERSE_DEBUG && !continue_flag && println("Infeasible node encountered.")

             # :^
            elseif op === 4
                REVERSE_DEBUG && println("power_rev[$k]      at k = $k -> $(setstorage[k])")
                continue_flag &= reverse_power!(k, children_arr, children_idx, numvalued, numberstorage,
                                                setstorage, x, lbd, ubd, sparsity, subgrad_tol, is_post)
                REVERSE_DEBUG && !continue_flag && println("Infeasible node encountered.")

            # :/
            elseif op === 5
                REVERSE_DEBUG && println("power_div[$k]      at k = $k -> $(setstorage[k])")
                continue_flag &= reverse_divide!(k, children_arr, children_idx, numvalued, numberstorage,
                                                 setstorage, x, lbd, ubd, sparsity, subgrad_tol, is_post)
                REVERSE_DEBUG && !continue_flag && println("Infeasible node encountered.")

            # ifelse
            elseif op === 6
                continue
            end

        # assumes that child is set-valued and thus parent is set-valued (since isnumber already checked)
        elseif nod.nodetype == JuMP._Derivatives.CALLUNIVAR
            op = nod.index
            REVERSE_DEBUG && println("fop_rev[$op][$k]       at k = $k -> $(setstorage[k])")
            if op <= JuMP._Derivatives.USER_UNIVAR_OPERATOR_ID_START
                arg_indx =  children_arr[adj.colptr[k]]
                continue_flag &= reverse_univariate!(k, op, arg_indx, setstorage, x, lbd, ubd, sparsity,
                                                     subgrad_tol, is_post)
            end
            REVERSE_DEBUG && !continue_flag && println("Infeasible node encountered.")
        end

        !continue_flag && break
    end

    return continue_flag
end

"""
$(FUNCTIONNAME)

A reverse_pass! on a `BufferedNonlinear` structure `d` intersects the existing value of the `d` with
constraint bounds then reverse propagates a set-valued operator (by default McCormick operator) along the
computational tape. The tapes are updated in place and boolean value is returned indicating whether the
reverse propagation yeilded a infeasible point (true = still feasible, false is proved infeasible).
"""
function reverse_pass!(evaluator::Evaluator, d::NonlinearExpression{V}) where V

    return reverse_pass_kernel!(d.nd, d.adj, evaluator.x, evaluator.lower_variable_bounds,
                                evaluator.upper_variable_bounds, d.grad_sparsity,
                                d.setstorage, evaluator.subgrad_tol, d.numberstorage,
                                d.isnumber, evaluator.reverse_subgrad_tighten)
end

function reverse_pass!(evaluator::Evaluator, d::BufferedNonlinearFunction{V}) where V
    d.last_past_reverse = true
    set_intersect_value!(d.expr, Interval(d.lower_bound, d.upper_bound))
    return reverse_pass!(evaluator, d.expr)
end
