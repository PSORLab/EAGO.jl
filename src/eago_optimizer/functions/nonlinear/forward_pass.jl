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
# src/eago_optimizer/evaluator/passes.jl
# Functions used to compute forward and reverse pass of nonlinear functions.
#############################################################################

"""
$(FUNCTIONNAME)

Post process set_value operator. By default, performs the affine interval cut on
a MC structure.
"""
function set_value_post(x_values::Vector{Float64}, val::MC{N,T}, lower_variable_bounds::Vector{Float64},
                        upper_variable_bounds::Vector{Float64}) where {N, T<:RelaxTag}

    lower = val.cv
    upper = val.cc
    lower_refinement = true
    upper_refinement = true

    for i = 1:N

        x_val = @inbounds x_values[i]
        cv_val = @inbounds val.cv_grad[i]
        cc_val = @inbounds val.cc_grad[i]

        lower_bound = @inbounds lower_variable_bounds[i]
        upper_bound = @inbounds upper_variable_bounds[i]

        if lower_refinement
            if cv_val > 0.0
                if isinf(lower_bound)
                    !upper_refinement && break
                    lower_refinement = false
                else
                    lower += cv_val*(lower_bound - x_val)
                end
            else
                if isinf(upper_bound)
                    !upper_refinement && break
                    lower_refinement = false
                else
                    lower += cv_val*(upper_bound - x_val)
                end
            end
        end

        if upper_refinement
            if cc_val > 0.0
                if isinf(lower_bound)
                    !lower_refinement && break
                    upper_refinement = false
                else
                    upper += cc_val*(upper_bound - x_val)
                end
            else
                if isinf(upper_bound)
                    !lower_refinement && break
                    upper_refinement = false
                else
                    upper += cc_val*(lower_bound - x_val)
                end
            end
        end
    end

    lower = lower_refinement ? max(lower, val.Intv.lo) : val.Intv.lo
    upper = upper_refinement ? min(upper, val.Intv.hi) : val.Intv.hi

    return MC{N,T}(val.cv, val.cc, Interval{Float64}(lower, upper), val.cv_grad, val.cc_grad, val.cnst)
end

"""
$(FUNCTIONNAME)

Intersects the new set valued operator with the prior and performing `set_value_post` if the flag is...
"""
function overwrite_or_intersect(xMC::MC{N,T}, past_xMC::MC{N,T}, x::Vector{Float64}, lbd::Vector{Float64},
                                ubd::Vector{Float64}, is_post::Bool, is_intersect::Bool)

    if is_post && is_intersect
        return set_value_post(x, xMC ∩ past_xMC, lbd, ubd)
    elseif is_post && !is_intersect
        return set_value_post(x, xMC, lbd, ubd)
    elseif !is_post && is_intersect
        return x ∩ past_xMC
    end
    return xMC
end

#=
macro get_binary_storage(children_idx, children_arr, numvalued, numberstorage, setstorage, type)
    esc(quote
        # get row indices
        idx1 = first($children_idx)
        idx2 = last($children_idx)

        # extract values for argument 1
        arg1_index = @inbounds $children_arr[idx1]
        arg1_is_number = @inbounds $numvalued[arg1_index]
        if arg1_is_number
            set1 = zero($type)
            num1 = @inbounds $numberstorage[arg1_index]
        else
            num1 = 0.0
            set1 = @inbounds $setstorage[arg1_index]
        end

        # extract values for argument 2
        arg2_index = @inbounds $children_arr[idx2]
        arg2_is_number = @inbounds $numvalued[arg2_index]
        if arg2_is_number
            num2 = @inbounds $numberstorage[arg2_index]
            set2 = zero(type)
        else
            set2 = @inbounds $setstorage[arg2_index]
            num2 = 0.0
        end
    end)
end
=#

function forward_plus_binary!()

    # get row indices
    idx1 = first(children_idx)
    idx2 = last(children_idx)

    # extract values for argument 1
    arg1_index = @inbounds children_arr[idx1]
    arg1_is_number = @inbounds numvalued[arg1_index]
    if arg1_is_number
        set1 = zero(MC{N,T})
        num1 = @inbounds numberstorage[arg1_index]
    else
        num1 = 0.0
        set1 = @inbounds setstorage[arg1_index]
    end

    # extract values for argument 2
    arg2_index = @inbounds children_arr[idx2]
    arg2_is_number = @inbounds numvalued[arg2_index]
    if arg2_is_number
        num2 = @inbounds numberstorage[arg2_index]
        set2 = zero(MC{N,T})
    else
        set2 = @inbounds setstorage[arg2_index]
        num2 = 0.0
    end

    output_is_number = arg1_is_number && arg2_is_number

    # a + b
    if output_is_number
        @inbounds numberstorage[k] = num1 + num2

    # x + b
    elseif !arg1_is_number && arg2_is_number
        outset = is_first_eval ? (set1 + num2) : plus_kernel(set1, num2, setstorage[k].Intv)

    # a + y
    elseif arg1_is_number && !arg2_is_number
        outset = is_first_eval ? (num1 + set2) : plus_kernel(num1, set2, setstorage[k].Intv)

    # x + y
    else
        outset = is_first_eval ? (set1 + set2) : plus_kernel(set1, set2, setstorage[k].Intv)

    end

    @inbounds numvalued[k] = output_is_number
    if !isnum
        inbounds setstorage[k] = overwrite_or_intersect(outset, setstorage[k], x_values, lbd, ubd, is_post, is_intersect)
    end

    return nothing
end

function forward_plus_narity!(k::Int64, children_idx::UnitRange{Int64}, children_arr::Vector{Int64},
                       numvalued::Vector{Bool}, numberstorage::Vector{Float64}, setstorage::Vector{MC{N,T}},
                       first_eval_flag::Bool) where {N, T<:RelaxTag}


    # get row indices
    idx = first(children_idx)

    # extract values for argument 1
    arg_index = @inbounds children_arr[idx]
    output_is_number = @inbounds numvalued[arg1_index]
    if output_is_number
       tmp_set = zero(MC{N,T})
       tmp_num = @inbounds numberstorage[arg_index]
    else
       tmp_num = 0.0
       tmp_set = @inbounds setstorage[arg_index]
    end

    output_is_number = true

    for idx = 2:length(children_idx)
        arg_index = @inbounds children_arr[idx]
        arg_is_number = @inbounds numvalued[arg_index]
        if arg_is_number
            tmp_num += @inbounds numberstorage[arg_index]
        else
            tmp_set += @inbounds setstorage[arg_index]
        end
        output_is_number &= arg_is_number
    end

    @inbounds numvalued[k] = output_is_number
    if isnum
        @inbounds numberstorage[k] = tmp_num
    else
        tmp_set += tmp_num
        setstorage[k] = overwrite_or_intersect(tmp_set, setstorage[k], x_values, lbd, ubd, is_post, is_intersect)
    end

    return nothing
end

function forward_multiply_binary!()
    # get row indices
    idx1 = first(children_idx)
    idx2 = last(children_idx)

    # extract values for argument 1
    arg1_index = @inbounds children_arr[idx1]
    arg1_is_number = @inbounds numvalued[arg1_index]
    if arg1_is_number
        set1 = zero(MC{N,T})
        num1 = @inbounds numberstorage[arg1_index]
    else
        num1 = 0.0
        set1 = @inbounds setstorage[arg1_index]
    end

    # extract values for argument 2
    arg2_index = @inbounds children_arr[idx2]
    arg2_is_number = @inbounds numvalued[arg2_index]
    if arg2_is_number
        num2 = @inbounds numberstorage[arg2_index]
        set2 = zero(MC{N,T})
    else
        set2 = @inbounds setstorage[arg2_index]
        num2 = 0.0
    end

    output_is_number = arg1_is_number && arg2_is_number

    # a * b
    if output_is_number
        @inbounds numberstorage[k] = num1 * num2

    # x * b
    elseif !arg1_is_number && arg2_is_number
        outset = is_first_eval ? (set1 * num2) : mult_kernel(set1, num2, setstorage[k].Intv)

    # a * y
    elseif arg1_is_number && !arg2_is_number
        outset = is_first_eval ? (num1 * set2) : mult_kernel(num1, set2, setstorage[k].Intv)

    # x * y
    else
        outset = is_first_eval ? (set1 * set2) : mult_kernel(set1, set2, setstorage[k].Intv)

    end

    @inbounds numvalued[k] = output_is_number
    if !isnum
        inbounds setstorage[k] = overwrite_or_intersect(outset, setstorage[k], x_values, lbd, ubd, is_post, is_intersect)
    end

    return nothing
end

function forward_multiply_narity!(k::Int64, x_values::Vector{Float64}, children_idx::UnitRange{Int64}, children_arr::Vector{Int64},
                           numvalued::Vector{Bool}, numberstorage::Vector{Float64}, setstorage::Vector{MC{N,T}},
                           current_node::NodeBB, subgrad_tighten::Bool, first_eval_flag::Bool) where {N, T<:RelaxTag}

    # get row indices
    idx = first(children_idx)

    # extract values for argument 1
    arg_index = @inbounds children_arr[idx]
    output_is_number = @inbounds numvalued[arg1_index]
    if output_is_number
        tmp_set = zero(MC{N,T})
        tmp_num = @inbounds numberstorage[arg_index]
    else
        tmp_num = 0.0
        tmp_set = @inbounds setstorage[arg_index]
    end

    output_is_number = true

    for idx = 2:length(children_idx)
        arg_index = @inbounds children_arr[idx]
        arg_is_number = @inbounds numvalued[arg_index]
        if arg_is_number
            tmp_num *= @inbounds numberstorage[arg_index]
        else
            tmp_set *= @inbounds setstorage[arg_index]
        end
        output_is_number &= arg_is_number
    end

    @inbounds numvalued[k] = output_is_number
    if isnum
        @inbounds numberstorage[k] = tmp_num
    else
       tmp_set *= tmp_num
       setstorage[k] = overwrite_or_intersect(tmp_set, setstorage[k], x_values, lbd, ubd, is_post, is_intersect)
    end

    return nothing
end

function forward_minus!(k::Int64, x_values::Vector{Float64}, children_idx::UnitRange{Int64},
                        children_arr::Vector{Int64}, numvalued::Vector{Bool}, numberstorage::Vector{Float64},
                        setstorage::Vector{MC{N,T}}, current_node::NodeBB, subgrad_tighten::Bool,
                        is_first_eval::Bool) where {N, T<:RelaxTag}

    # get row indices
    idx1 = first(children_idx)
    idx2 = last(children_idx)

    # extract values for argument 1
    arg1_index = @inbounds children_arr[idx1]
    arg1_is_number = @inbounds numvalued[arg1_index]
    if arg1_is_number
        set1 = zero(MC{N,T})
        num1 = @inbounds numberstorage[arg1_index]
    else
        num1 = 0.0
        set1 = @inbounds setstorage[arg1_index]
    end

    # extract values for argument 2
    arg2_index = @inbounds children_arr[idx2]
    arg2_is_number = @inbounds numvalued[arg2_index]
    if arg2_is_number
        num2 = @inbounds numberstorage[arg2_index]
        set2 = zero(MC{N,T})
    else
        set2 = @inbounds setstorage[arg2_index]
        num2 = 0.0
    end

    output_is_number = arg1_is_number && arg2_is_number

    # a - b
    if output_is_number
        @inbounds numberstorage[k] = num1 - num2

    # x - b
    elseif !arg1_is_number && arg2_is_number
        outset = is_first_eval ? (set1 - num2) : minus_kernel(set1, num2, setstorage[k].Intv)

    # a - y
    elseif arg1_is_number && !arg2_is_number
        outset = is_first_eval ? (num1 - set2) : minus_kernel(num1, set2, setstorage[k].Intv)

    # x - y
    else
        outset = is_first_eval ? (set1 - set2) : minus_kernel(set1, set2, setstorage[k].Intv)

    end

    @inbounds numvalued[k] = output_is_number
    if ~isnum
        inbounds setstorage[k] = overwrite_or_intersect(outset, setstorage[k], x_values, lbd, ubd, is_post, is_intersect)
    end

    return nothing
end

function forward_power!(k::Int64, x_values::Vector{Float64}, children_idx::UnitRange{Int64},
                        children_arr::Vector{Int64}, numvalued::Vector{Bool}, numberstorage::Vector{Float64},
                        setstorage::Vector{MC{N,T}}, current_node::NodeBB, subgrad_tighten::Bool,
                        is_first_eval::Bool,
                        ctx::GuardCtx) where {N, T<:RelaxTag}

    # get row indices
    idx1 = first(children_idx)
    idx2 = last(children_idx)

    # extract values for argument 1
    arg1_index = @inbounds children_arr[idx1]
    arg1_is_number = @inbounds numvalued[arg1_index]
    if arg1_is_number
        set1 = zero(MC{N,T})
        num1 = @inbounds numberstorage[arg1_index]
    else
        num1 = 0.0
        set1 = @inbounds setstorage[arg1_index]
    end

    # extract values for argument 2
    arg2_index = @inbounds children_arr[idx2]
    arg2_is_number = @inbounds numvalued[arg2_index]
    if arg2_is_number
        num2 = @inbounds numberstorage[arg2_index]
        set2 = zero(MC{N,T})
    else
        set2 = @inbounds setstorage[arg2_index]
        num2 = 0.0
    end

    # is output a number (by closure of the reals)?
    output_is_number = arg1_is_number && arg2_is_number
    @inbounds numvalued[k] = output_is_number

    # x^1 = x
    if num2 === 1.0
        if arg1_is_number
            @inbounds numberstorage[k] = num1
        else
            @inbounds setstorage[k] = set1
        end
        return nothing

    # x^0 = 1
    elseif num2 === 0.0
        if arg1_is_number
            @inbounds numberstorage[k] = 1.0
        else
            @inbounds setstorage[k] = zero(MC{N,T})
        end
        return nothing

    else
        # a^b
        if arg1_is_number && arg2_is_number
            @inbounds numberstorage[k] = num1^num2

        # x^b
        elseif !arg1_is_number && arg2_is_number
            outset = is_first_eval ? pow(set1, num2) : ^(set1, num2, setstorage[k].Intv)

        # a^y
        elseif arg1_is_number  && !arg2_is_number
            outset = is_first_eval ? overdub(ctx, pow, num1, set2) : overdub(ctx, ^, num1, set2, setstorage[k].Intv)

        # x^y
        elseif !arg1_is_number && !arg2_is_number
            outset = is_first_eval ? overdub(ctx, pow, set1, set2) : overdub(ctx, ^, set1, set2, setstorage[k].Intv)

        end
    end

    if !output_is_number
        setstorage[k] = overwrite_or_intersect(outset, setstorage[k], x_values, lbd, ubd, is_post, is_intersect)
    end

    return nothing
end

function forward_divide!(k::Int64, x_values::Vector{Float64}, children_idx::UnitRange{Int64},
                         children_arr::Vector{Int64}, numvalued::Vector{Bool}, numberstorage::Vector{Float64},
                         setstorage::Vector{MC{N,T}},
                         current_node::NodeBB,
                         subgrad_tighten::Bool,
                         first_eval_flag::Bool,
                         ctx::GuardCtx) where {N,T<:RelaxTag}

    # get row indices
    idx1 = first(children_idx)
    idx2 = last(children_idx)

    # extract values for argument 1
    arg1_index = @inbounds children_arr[idx1]
    arg1_is_number = @inbounds numvalued[arg1_index]
    if arg1_is_number
        set1 = zero(MC{N,T})
        num1 = @inbounds numberstorage[arg1_index]
    else
        num1 = 0.0
        set1 = @inbounds setstorage[arg1_index]
    end

    # extract values for argument 2
    arg2_index = @inbounds children_arr[idx2]
    arg2_is_number = @inbounds numvalued[arg2_index]
    if arg2_is_number
        num2 = @inbounds numberstorage[arg2_index]
        set2 = zero(MC{N,T})
    else
        set2 = @inbounds setstorage[arg2_index]
        num2 = 0.0
    end

    # is output a number (by closure of the reals)?
    output_is_number = arg1_is_number && arg2_is_number
    @inbounds numvalued[k] = output_is_number

    # a/b
    if output_is_number
        @inbounds numberstorage[k] = num1/num2

    # x/b
    elseif !arg1_is_number && arg2_is_number
        outset = is_first_eval ? set1/num2 : div_kernel(set1, num2, setstorage[k].Intv)

    # a/y
    elseif arg1_is_number && !arg2_is_number
        outset = is_first_eval ? num1/set2 : div_kernel(num1, set2, setstorage[k].Intv)

    # x/y
    else
        outset = is_first_eval ? set1/set2 : div_kernel(set1, set2, setstorage[k].Intv)

    end

    @inbounds setstorage[k] = overwrite_or_intersect(outset, setstorage[k], x_values, lbd, ubd, is_post, is_intersect)

    return nothing
end

function forward_user_multivariate!()
        op_sym = id_to_operator[op]
        evaluator = user_operators.multivariate_operator_evaluator[op - JuMP._Derivatives.USER_OPERATOR_ID_START+1]
        f_input = view(user_input_buffer, 1:n_children)
        fnum_input = view(flt_user_input_buffer, 1:n_children)
        r = 1
        isnum = true
        for c_idx in children_idx
            @inbounds ix = children_arr[c_idx]
            @inbounds chdset = numvalued[ix]
            isnum &= chdset
            if chdset
                if isnum
                    @inbounds fnum_input[r] = numberstorage[ix]
                end
                @inbounds f_input[r] = numberstorage[ix]
            else
                @inbounds f_input[r] = setstorage[ix]
            end
            r += 1
        end
        if isnum
            numberstorage[k] = MOI.eval_objective(evaluator, fnum_input)
        else
            fval = MOI.eval_objective(evaluator, f_input)
            #fval = Cassette.overdub(ctx, MOI.eval_objective, evaluator, f_input)
            setstorage[k] = overwrite_or_intersect(fval, setstorage[k], x, lbd, ubd, is_post, is_intersect)
        end
        numvalued[k] = isnum
    else
        error("Unsupported operation $(operators[op])")
    end

    return nothing
end

function forward_univariate_number!()
    return nothing
end

function forward_univariate_tiepnt_1!()
    @inbounds tpdict_tuple = tpdict[k]
    tindx1 = tpdict_tuple[1]
    tindx2 = tpdict_tuple[2]
    @inbounds tp1 = tp1storage[tindx1]
    @inbounds tp2 = tp2storage[tindx2]
    fval_mc, tp1, tp2 = Cassette.overdub(ctx, single_tp_set, op, child_val_mc, setstorage[k], tp1, tp2, first_eval_flag)
    @inbounds tp1storage[tindx] = tp1
    @inbounds tp1storage[tindx] = tp2
    setstorage[k] = overwrite_or_intersect(fval_mc, setstorage[k], x_values, lbd, ubd, is_post, is_intersect)
    return nothing
end

function forward_univariate_tiepnt_2!()
    @inbounds tpdict_tuple = tpdict[k]
    tindx1 = tpdict_tuple[1]
    tindx2 = tpdict_tuple[2]
    tindx3 = tpdict_tuple[3]
    tindx4 = tpdict_tuple[4]
    @inbounds tp1 = tp1storage[tindx1]
    @inbounds tp2 = tp2storage[tindx2]
    @inbounds tp3 = tp1storage[tindx3]
    @inbounds tp4 = tp2storage[tindx4]
    fval_mc, tp1, tp2, tp3, tp4 = Cassette.overdub(ctx, double_tp_set, op, child_val_mc, setstorage[k], tp1, tp2, tp3, tp4, first_eval_flag)
    @inbounds tp1storage[tindx1] = tp1
    @inbounds tp2storage[tindx2] = tp2
    @inbounds tp3storage[tindx1] = tp3
    @inbounds tp4storage[tindx2] = tp4
    setstorage[k] = overwrite_or_intersect(fval_mc, setstorage[k], x_values, lbd, ubd, is_post, is_intersect)
    return nothing
end

function forward_univariate_user!()
    userop = op - JuMP._Derivatives.USER_UNIVAR_OPERATOR_ID_START + 1
    @inbounds f = user_operators.univariate_operator_f[userop]
    if arg_is_number
        fval_num = f(child_val_num)
        @inbounds numberstorage[k] = fval_num
    else
        fval_mc = f(child_val_mc) #Cassette.overdub(ctx, f, child_val_mc)
        setstorage[k] = overwrite_or_intersect(fval_mc, setstorage[k], x_values, lbd, ubd, is_post, is_intersect)
    end
    return nothing
end

function forward_univariate_other!()
    if chdset
        fval_num = eval_univariate_set(op, child_val_num)
        @inbounds numberstorage[k] = fval_num
    else
        #eval_univariate_set(op, child_val_mc)
        fval_mc = Cassette.overdub(ctx, eval_univariate_set, op, child_val_mc)
        setstorage[k] = overwrite_or_intersect(fval_mc, setstorage[k], x_values, lbd, ubd, is_post, is_intersect)
    end
    return nothing
end

const id_to_operator = Dict(value => key for (key, value) in JuMP.univariate_operator_to_id)
function forward_pass_kernel!(setstorage::Vector{MC{N,T}}, numberstorage::Vector{Float64}, numvalued::Vector{Bool},
                       nd::Vector{JuMP.NodeData}, adj::SparseMatrixCSC{Bool,Int64},
                       current_node::NodeBB, x_values::Vector{Float64},
                       subexpr_values_flt::Vector{Float64}, subexpr_values_set::Vector{MC{N,T}},
                       subexpression_isnum::Vector{Bool}, user_input_buffer::Vector{MC{N,T}}, flt_user_input_buffer::Vector{Float64},
                       subgrad_tighten::Bool,
                       tpdict::Dict{Int64,Tuple{Int64,Int64,Int64,Int64}}, tp1storage::Vector{Float64},
                       tp2storage::Vector{Float64}, tp3storage::Vector{Float64}, tp4storage::Vector{Float64},
                       first_eval_flag::Bool, user_operators::JuMP._Derivatives.UserOperatorRegistry,
                       ctx::GuardCtx) where {N,T<:RelaxTag}

    children_arr = rowvals(adj)

    for k = length(nd):-1:1
        @inbounds nod = nd[k]
        op = nod.index

        if nod.nodetype == JuMP._Derivatives.VALUE
            numvalued[k] = true

        elseif nod.nodetype == JuMP._Derivatives.PARAMETER
            numvalued[k] = true

        elseif nod.nodetype == JuMP._Derivatives.VARIABLE
            xMC = MC{N,T}(x_values[op], Interval{Float64}(lvb[op], uvb[op]), Val(N))

            if first_eval_flag
                @inbounds setstorage[k] = xMC
            else
                @inbounds setstorage[k] = xMC ∩ setstorage[k]
            end
            numvalued[k] = false

        elseif nod.nodetype == JuMP._Derivatives.SUBEXPRESSION                          # DONE
            @inbounds isnum = subexpression_isnum[op]
            if isnum
                @inbounds numberstorage[k] = subexpr_values_flt[op]
            else
                @inbounds setstorage[k] = subexpr_values_set[op] #∩ setstorage[k]
            end
            numvalued[k] = isnum

        elseif nod.nodetype == JuMP._Derivatives.CALL

            @inbounds children_idx = nzrange(adj,k)
            n_children = length(children_idx)

            # :+ with arity two or greater
            if op === 1 # :+
                n = length(children_idx)
                if n === 2
                    forward_plus_binary!(k, children_idx, children_arr, numvalued, numberstorage, setstorage, first_eval_flag)
                else
                    forward_plus_narity!(k, children_idx, children_arr, numvalued, numberstorage, setstorage, first_eval_flag)
                end

            # :- with arity two
            elseif op === 2
                forward_minus!(k, x_values, children_idx, children_arr, numvalued, numberstorage,
                              setstorage, current_node, subgrad_tighten, first_eval_flag)

            # :* with arity two or greater
            elseif op === 3
                n = length(children_idx)
                if n === 2
                    forward_multiply_binary!(k, x_values, children_idx, children_arr, numvalued, numberstorage,
                                             setstorage, current_node, subgrad_tighten, first_eval_flag)
                else
                    forward_multiply_narity!(k, x_values, children_idx, children_arr, numvalued, numberstorage,
                                             setstorage, current_node, subgrad_tighten, first_eval_flag)
                end

            # :^
            elseif op === 4
                forward_power!(k, x_values, children_idx, children_arr, numvalued, numberstorage,
                               setstorage, current_node, subgrad_tighten, first_eval_flag, ctx)

            # :/
            elseif op === 5
                forward_divide!(k, x_values, children_idx, children_arr, numvalued, numberstorage,
                               setstorage, current_node, subgrad_tighten, first_eval_flag, ctx)

            elseif op >= JuMP._Derivatives.USER_OPERATOR_ID_START
                forward_user_multivariate!()

        elseif nod.nodetype == JuMP._Derivatives.CALLUNIVAR

            # checks to see if operator is a number
            @inbounds arg_idx = children_arr[adj.colptr[k]]
            arg_is_number = @inbounds numvalued[arg_idx]
            @inbounds numvalued[k] = arg_is_number

            # performs univariate operators on number valued inputs
            if arg_is_number
                forward_univariate_number!()

            # performs set valued operators that require a single tiepoint calculation
            elseif single_tp(op)
                forward_univariate_tiepnt_1!()

            # performs set valued operators that require two tiepoint calculations
            elseif double_tp(op)
                forward_univariate_tiepnt_2!()

            # performs set valued operators on user-defined univariate functions
            elseif op >= JuMP._Derivatives.USER_UNIVAR_OPERATOR_ID_START
                forward_univariate_user!()

            else
                forward_univariate_other!()

            end

        else
            error("Unrecognized node type $(nod.nodetype).")

        end
    end

    return nothing
end
