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

function forward_plus!(k::Int64, children_idx::UnitRange{Int64}, children_arr::Vector{Int64},
                       numvalued::Vector{Bool}, numberstorage::Vector{Float64}, setstorage::Vector{MC{N,T}},
                       first_eval_flag::Bool) where {N, T<:RelaxTag}

    #println("k = $k, numvalued[$k] = $(numvalued[k])")
    tmp_num = 0.0
    tmp_mc = zero(MC{N,T})
    isnum = true
    chdset = true
    n = length(children_idx)
    if n == 2
        @inbounds c_idx_1 = children_idx[1]
        @inbounds c_idx_2 = children_idx[2]
        @inbounds ix1 = children_arr[c_idx_1]
        @inbounds ix2 = children_arr[c_idx_2]
        @inbounds chdset1 = numvalued[ix1]
        @inbounds chdset2 = numvalued[ix2]
        if first_eval_flag
            if chdset1
                tmp_num += numberstorage[ix1]
            else
                tmp_mc += setstorage[ix1]
            end
            if chdset2
                tmp_num += numberstorage[ix2]
            else
                tmp_mc += setstorage[ix2]
            end
            isnum = (chdset1 & chdset2)
            @inbounds numvalued[k] = isnum
        else
            if ~chdset1 && ~chdset2
                tmp_mc = McCormick.plus_kernel(setstorage[ix1], setstorage[ix2], setstorage[k].Intv)
            elseif chdset1 && ~chdset2
                tmp_mc = McCormick.plus_kernel(numberstorage[ix1], setstorage[ix2], setstorage[k].Intv)
            elseif chdset2 && ~chdset1
                tmp_mc = McCormick.plus_kernel(setstorage[ix1], numberstorage[ix2], setstorage[k].Intv)
            else
                @inbounds tmp_num = numberstorage[k]
            end
        end
    else
        for c_idx in children_idx
            @inbounds ix = children_arr[c_idx]
            @inbounds chdset = numvalued[ix]
            if chdset
                @inbounds tmp_num += numberstorage[ix]
            else
                @inbounds tmp_mc += setstorage[ix]
            end
            isnum &= chdset
        end
        @inbounds numvalued[k] = isnum
    end
    if isnum
        @inbounds numberstorage[k] = tmp_num
    else
        setstorage[k] = overwrite_or_intersect(tmp_mc_1, setstorage[k], x_values, lbd, ubd, is_post, is_intersect)
    end

    return
end

function forward_minus!(k::Int64, x_values::Vector{Float64}, ix1::Int64, ix2::Int64, numvalued::Vector{Bool},
                        numberstorage::Vector{Float64}, setstorage::Vector{MC{N,T}},
                        current_node::NodeBB, subgrad_tighten::Bool, first_eval_flag::Bool) where {N, T<:RelaxTag}

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

function forward_multiply!(k::Int64, x_values::Vector{Float64}, children_idx::UnitRange{Int64}, children_arr::Vector{Int64},
                           numvalued::Vector{Bool}, numberstorage::Vector{Float64}, setstorage::Vector{MC{N,T}},
                           current_node::NodeBB, subgrad_tighten::Bool, first_eval_flag::Bool) where {N, T<:RelaxTag}
    #println("forward multiply kernel")
    tmp_num_1 = 1.0
    tmp_mc_1 = one(MC{N,T})
    isnum = true
    chdset = true
    is_first = true
    for c_idx in children_idx
        cix = @inbounds children_arr[c_idx]
        chdset = @inbounds numvalued[cix]
        isnum &= chdset
        if chdset
            tmp_num_1 *= @inbounds numberstorage[cix]
        else
            if is_first
                tmp_mc_1 = @inbounds setstorage[cix]
                is_first = false
            else
                tmp_mc_1 = tmp_mc_1*(@inbounds setstorage[cix])
            end
        end
    end
    if isnum
        numberstorage[k] = tmp_num_1
    else
        tmp_mc_1 *= tmp_num_1
        setstorage[k] = overwrite_or_intersect(tmp_mc_1, setstorage[k], x_values, lbd, ubd, is_post, is_intersect)
    end
    numvalued[k] = isnum

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

    # x^1 = x
    if num2 === 1.0
        if arg1_is_number
            @inbounds numberstorage[k] = num1
        else
            @inbounds setstorage[k] = set1
        end

    # x^0 = 1
    elseif num2 === 0.0
        if arg1_is_number
            @inbounds numberstorage[k] = 1.0
        else
            @inbounds setstorage[k] = zero(MC{N,T})
        end

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

    @inbounds numvalued[k] = output_is_number
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
    tmp_num_1 = 0.0
    tmp_num_2 = 0.0
    tmp_mc_1 = zero(MC{N,T})
    tmp_mc_2 = zero(MC{N,T})
    idx1 = first(children_idx)
    idx2 = last(children_idx)
    @inbounds ix1 = children_arr[idx1]
    @inbounds ix2 = children_arr[idx2]
    @inbounds chdset1 = numvalued[ix1]
    @inbounds chdset2 = numvalued[ix2]
    if chdset1
        @inbounds tmp_num_1 = numberstorage[ix1]
    else
        @inbounds tmp_mc_1 = setstorage[ix1]
    end
    if chdset2
        @inbounds tmp_num_2 = numberstorage[ix2]
    else
        @inbounds tmp_mc_2 = setstorage[ix2]
    end
    if chdset1 && chdset2
        @inbounds numberstorage[k] = tmp_num_1/tmp_num_2
    else
        if first_eval_flag
            if !chdset1 && chdset2
                tmp_mc_1 = tmp_mc_1/tmp_num_2
            elseif chdset1 && !chdset2
                tmp_mc_1 = tmp_num_1/tmp_mc_2
            elseif !chdset1 && !chdset2
                tmp_mc_1 = tmp_mc_1/tmp_mc_2 #TODO: Overdub this with guard context
            end
        else
            @inbounds intv_k = setstorage[k].Intv
            if !chdset1 && chdset2
                tmp_mc_1 = McCormick.div_kernel(tmp_mc_1, tmp_num_2, intv_k)
            elseif chdset1 && !chdset2
                tmp_mc_1 = McCormick.div_kernel(tmp_num_1, tmp_mc_2, intv_k)
            elseif !chdset1 && !chdset2
                tmp_mc_1 = McCormick.div_kernel(tmp_mc_1, tmp_mc_2, intv_k) #TODO: Overdub this with guard context
            end
        end
        setstorage[k] = overwrite_or_intersect(tmp_mc_1, setstorage[k], x_values, lbd, ubd, is_post, is_intersect)
    end
    @inbounds numvalued[k] = chdset1 & chdset2
    return
end

const id_to_operator = Dict(value => key for (key, value) in JuMP.univariate_operator_to_id)
function forward_eval(setstorage::Vector{MC{N,T}}, numberstorage::Vector{Float64}, numvalued::Vector{Bool},
                      nd::Vector{JuMP.NodeData}, adj::SparseMatrixCSC{Bool,Int64},
                      current_node::NodeBB, x_values::Vector{Float64},
                      subexpr_values_flt::Vector{Float64}, subexpr_values_set::Vector{MC{N,T}},
                      subexpression_isnum::Vector{Bool}, user_input_buffer::Vector{MC{N,T}}, flt_user_input_buffer::Vector{Float64},
                      subgrad_tighten::Bool,
                      tpdict::Dict{Int64,Tuple{Int64,Int64,Int64,Int64}}, tp1storage::Vector{Float64},
                      tp2storage::Vector{Float64}, tp3storage::Vector{Float64}, tp4storage::Vector{Float64},
                      first_eval_flag::Bool, user_operators::JuMP._Derivatives.UserOperatorRegistry,
                      seeds::Vector{SVector{N,Float64}},
                      ctx::GuardCtx) where {N,T<:RelaxTag}

    @assert length(numberstorage) >= length(nd)
    @assert length(setstorage) >= length(nd)
    @assert length(numvalued) >= length(nd)

    set_value_sto = zero(MC{N,T})
    children_arr = rowvals(adj)
    n = length(x_values)

    tmp_num_1 = 0.0
    tmp_num_2 = 0.0
    tmp_mc_1 = zero(MC{N,T})
    tmp_mc_2 = zero(MC{N,T})
    for k = length(nd):-1:1
        @inbounds nod = nd[k]
        op = nod.index
        if nod.nodetype == JuMP._Derivatives.VALUE
            numvalued[k] = true
        elseif nod.nodetype == JuMP._Derivatives.PARAMETER
            numvalued[k] = true
        elseif nod.nodetype == JuMP._Derivatives.VARIABLE
            seed = seeds[op]
            xMC = MC{N,T}(x_values[op], x_values[op],
                        Interval{Float64}(current_node.lower_variable_bounds[op],
                                     current_node.upper_variable_bounds[op]),
                                     seed, seed, false)

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
            if op == 1 # :+
                forward_plus!(k, children_idx, children_arr, numvalued, numberstorage, setstorage, first_eval_flag)::Nothing
            elseif op == 2 # :-
                #println("FORWARD MINUS IN")
                @assert n_children == 2
                child1 = first(children_idx)
                @inbounds ix1 = children_arr[child1]
                @inbounds ix2 = children_arr[child1+1]
                forward_minus!(k, x_values, ix1, ix2, numvalued, numberstorage,
                               setstorage, current_node, subgrad_tighten, first_eval_flag)::Nothing
                #println("FORWARD MINUS END[$k] = $(setstorage[k])")
            elseif op == 3 # :*
                #println("FORWARD MULTIPLY IN")
                forward_multiply!(k, x_values, children_idx, children_arr, numvalued, numberstorage,
                                  setstorage, current_node, subgrad_tighten, first_eval_flag)::Nothing
                #println("FORWARD MULTIPLY END[$k] = $(setstorage[k])")
            elseif op == 4 # :^                                                      # DONE
                @assert n_children == 2
                forward_power!(k, x_values, children_idx, children_arr, numvalued, numberstorage,
                               setstorage, current_node, subgrad_tighten, first_eval_flag, ctx)::Nothing
            elseif op == 5 # :/                                                      # DONE
                @assert n_children == 2
                forward_divide!(k, x_values, children_idx, children_arr, numvalued, numberstorage,
                               setstorage, current_node, subgrad_tighten, first_eval_flag, ctx)::Nothing
            elseif op >= JuMP._Derivatives.USER_OPERATOR_ID_START
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
        elseif nod.nodetype == JuMP._Derivatives.CALLUNIVAR
            @inbounds child_idx = children_arr[adj.colptr[k]]
            @inbounds chdset = numvalued[child_idx]
            if chdset
                @inbounds child_val_num = numberstorage[child_idx]
            else
                @inbounds child_val_mc = setstorage[child_idx]
            end
            if op >= JuMP._Derivatives.USER_UNIVAR_OPERATOR_ID_START
                userop = op - JuMP._Derivatives.USER_UNIVAR_OPERATOR_ID_START + 1
                @inbounds f = user_operators.univariate_operator_f[userop]
                if chdset
                    fval_num = f(child_val_num)
                    @inbounds numberstorage[k] = fval_num
                else
                    fval_mc = f(child_val_mc) #Cassette.overdub(ctx, f, child_val_mc)
                    setstorage[k] = overwrite_or_intersect(fval_mc, setstorage[k], x_values, lbd, ubd, is_post, is_intersect)
                end
            elseif single_tp(op) && chdset
                @inbounds tpdict_tuple = tpdict[k]
                tindx1 = tpdict_tuple[1]
                tindx2 = tpdict_tuple[2]
                @inbounds tp1 = tp1storage[tindx1]
                @inbounds tp2 = tp2storage[tindx2]
                fval_mc, tp1, tp2 = Cassette.overdub(ctx, single_tp_set, op, child_val_mc, setstorage[k], tp1, tp2, first_eval_flag)
                @inbounds tp1storage[tindx] = tp1
                @inbounds tp1storage[tindx] = tp2
                setstorage[k] = overwrite_or_intersect(fval_mc, setstorage[k], x_values, lbd, ubd, is_post, is_intersect)
            elseif double_tp(op) && chdset
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
            else
                if chdset
                    fval_num = eval_univariate_set(op, child_val_num)
                    @inbounds numberstorage[k] = fval_num
                else
                    #eval_univariate_set(op, child_val_mc)
                    fval_mc = Cassette.overdub(ctx, eval_univariate_set, op, child_val_mc)
                    setstorage[k] = overwrite_or_intersect(fval_mc, setstorage[k], x_values, lbd, ubd, is_post, is_intersect)
                end
            end
            numvalued[k] = chdset
        else
            error("Unrecognized node type $(nod.nodetype).")
        end
    end

    return
end

# maximum number to perform reverse operation on associative term by summing
# and evaluating pairs remaining terms not reversed
const MAX_ASSOCIATIVE_REVERSE = 4

"""
$(TYPEDSIGNATURES)
"""
function reverse_eval(setstorage::Vector{T}, numberstorage, numvalued, subexpression_isnum,
                      subexpr_values_set, nd::Vector{JuMP.NodeData}, adj, x_values, current_node::NodeBB, subgrad_tighten::Bool) where T

    @assert length(setstorage) >= length(nd)
    @assert length(numberstorage) >= length(nd)
    @assert length(numvalued) >= length(nd)

    children_arr = rowvals(adj)
    N = length(x_values)

    tmp_hold = zero(T)
    continue_flag = true

    for k = 1:length(nd)
        @inbounds nod = nd[k]
        if (nod.nodetype == JuMP._Derivatives.VALUE ||
            nod.nodetype == JuMP._Derivatives.LOGIC ||
            nod.nodetype == JuMP._Derivatives.COMPARISON ||
            nod.nodetype == JuMP._Derivatives.PARAMETER ||
            nod.nodetype == JuMP._Derivatives.EXTRA )
            continue                  # DONE
        elseif nod.nodetype == JuMP._Derivatives.VARIABLE
            op = nod.index
            @inbounds current_node.lower_variable_bounds[op] = setstorage[k].Intv.lo
            @inbounds current_node.upper_variable_bounds[op] = setstorage[k].Intv.hi               # DONE
        elseif nod.nodetype == JuMP._Derivatives.SUBEXPRESSION
            #=
            @inbounds isnum = subexpression_isnum[nod.index]
            if ~isnum
                @inbounds subexpr_values_set[nod.index] = setstorage[k]
            end          # DONE
            =#
        elseif numvalued[k]
            continue                                             # DONE
        elseif nod.nodetype == JuMP._Derivatives.CALL
            op = nod.index
            #println("  reverse pass op = $op at k = $k")
            parent_index = nod.parent
            @inbounds children_idx = nzrange(adj,k)
            @inbounds parent_value = setstorage[k]
            n_children = length(children_idx)
            if (op >= JuMP._Derivatives.USER_OPERATOR_ID_START)
                # SKIPS USER DEFINE OPERATORS NOT BRIDGED INTO JuMP Tree
                # Representation
                continue
            elseif op == 1 # :+

                #println("+")
                lenx = length(children_idx)
                count = 0
                child_arr_indx = children_arr[children_idx]
                for c_idx in child_arr_indx
                    tmp_sum = 0.0
                    @inbounds inner_chdset = numvalued[c_idx]
                    if ~inner_chdset
                        if count < MAX_ASSOCIATIVE_REVERSE
                            for cin_idx in child_arr_indx
                                if cin_idx != c_idx
                                    if @inbounds numvalued[cin_idx]
                                        tmp_sum += @inbounds numberstorage[cin_idx]
                                    else
                                        tmp_sum += @inbounds setstorage[cin_idx]
                                    end
                                end
                            end
                            @inbounds tmp_hold = setstorage[c_idx]
                            #println(" ")
                            #println(" pre + ")
                            #println("parent_value: $(parent_value)")
                            #println("tmp_hold: $(tmp_hold)")
                            #println("tmp_sum: $(tmp_sum)")
                            pnew, xhold, xsum = plus_rev(parent_value, tmp_hold, tmp_sum)
                            #println("post + ")
                            #println("pnew: $pnew")
                            #println("xhold: $xhold")
                            #println("xsum: $xsum")
                            #println(" ")
                            if isempty(pnew) || isempty(xhold) || isempty(xsum)
                                continue_flag = false
                                break
                            end
                            if isnan(pnew)
                                pnew = interval_MC(parent_value)
                            end
                            if isnan(xhold)
                                pnew = interval_MC(tmp_hold)
                            end
                            setstorage[k] = set_value_post(x_values, pnew, current_node, subgrad_tighten)
                            setstorage[c_idx] = set_value_post(x_values, xhold, current_node, subgrad_tighten)
                        else
                            break
                        end
                    end
                    count += 1
                end
                !continue_flag && break

            elseif op == 2 # :-

                #println("-")
                child1 = first(children_idx)
                child2 = last(children_idx)
                @assert n_children == 2
                @inbounds ix1 = children_arr[child1]
                @inbounds ix2 = children_arr[child2]
                @inbounds chdset1 = numvalued[ix1]
                @inbounds chdset2 = numvalued[ix2]
                #println(" ")
                #println(" pre + ")
                #println("parent_value: $(parent_value)")
                if chdset1 && !chdset2
                    pnew, xnew, ynew = minus_rev(parent_value, numberstorage[ix1], setstorage[ix2])
                    #println("arg1: $(numberstorage[ix1])")
                    #println("arg2: $(setstorage[ix2])")
                    #println("post + ")
                    #println("pnew: $pnew")
                    #println("xhold: $xnew")
                    #println("xsum: $ynew")
                elseif !chdset1 && chdset2
                    pnew, xnew, ynew = minus_rev(parent_value, setstorage[ix1], numberstorage[ix2])
                    #println("arg1: $(setstorage[ix1])")
                    #println("arg2: $(numberstorage[ix2])")
                    #println("post + ")
                    #println("pnew: $pnew")
                    #println("xhold: $xnew")
                    #println("xsum: $ynew")
                elseif chdset1 && chdset2
                    pnew, xnew, ynew = minus_rev(parent_value, setstorage[ix1], setstorage[ix2])
                    #println("arg1: $(setstorage[ix1])")
                    #println("arg2: $(setstorage[ix2])")
                    #println("post + ")
                    #println("pnew: $pnew")
                    #println("xhold: $xnew")
                    #println("xsum: $ynew")
                end
                #println(" ")
                if !(!chdset1 && !chdset2)
                    flagp = isempty(pnew)
                    flagx = isempty(xnew)
                    flagy = isempty(ynew)
                    if (flagp || flagx || flagy)
                        continue_flag = false
                        break
                    end
                    if isnan(pnew)
                        pnew = interval_MC(parent_value)
                    end
                    setstorage[k] = pnew
                    if ~chdset1
                        if isnan(xnew)
                            xnew = interval_MC(setstorage[ix1])
                        end
                        #println("a1 pnew, xnew, ynew = $pnew, $xnew, $ynew)"
                        setstorage[ix1] = set_value_post(x_values, xnew, current_node, subgrad_tighten)
                    end
                    if  ~chdset2
                        if isnan(ynew)
                            ynew = interval_MC(setstorage[ix2])
                        end
                        #println("a2 pnew, xnew, ynew = $pnew, $xnew, $ynew)"
                        setstorage[ix2] = set_value_post(x_values, ynew, current_node, subgrad_tighten)
                    end
                end
                #println("setstorage[$k] = $(setstorage[k]), setstorage[$ix1] = $(setstorage[ix1]), setstorage[$ix2] = $(setstorage[ix2])")

            elseif (op == 3) # :*
                #println(" ")
                #println("start *")

                tmp_mlt = 1.0
                chdset = true
                count = 0
                child_arr_indx = children_arr[children_idx]
                for c_idx in child_arr_indx
                    if count < MAX_ASSOCIATIVE_REVERSE
                        if ~numvalued[c_idx]
                            tmp_mlt = 1.0
                            for cin_idx in child_arr_indx
                                if cin_idx != c_idx
                                    @inbounds chdset = numvalued[cin_idx]
                                    #println("numvalued[cin_idx] = $(numvalued[cin_idx])")
                                    #println("numberstorage[cin_idx] = $(numberstorage[cin_idx])")
                                    #println("setstorage[cin_idx] = $(setstorage[cin_idx])")
                                    if chdset
                                        @inbounds tmp_mlt *= numberstorage[cin_idx]
                                    else
                                        @inbounds tmp_mlt *= setstorage[cin_idx]
                                    end
                                end
                            end
                            @inbounds chdset = numvalued[c_idx]
                            #println("numvalued[c_idx] = $(numvalued[c_idx])")
                            if chdset
                                #println(" ")
                                #println(" input 1 ")
                                #println("parent_value: $(parent_value)")
                                #println("numberstorage[c_idx]: $(numberstorage[c_idx])")
                                #println("setstorage[c_idx]: $(setstorage[c_idx])")
                                #println("tmp_mlt: $(tmp_mlt)")
                                @inbounds pnew, xhold, xprd = mul_rev(parent_value, numberstorage[c_idx], tmp_mlt)
                                #println(" output ")
                                #println(" pnew = $pnew")
                                #println(" xhold = $xhold")
                                #println(" xprd = $xprd")
                            else
                                #println(" ")
                                #println(" input 2 ")
                                #println("parent_value: $(parent_value)")
                                #println("numberstorage[c_idx]: $(numberstorage[c_idx])")
                                #println("setstorage[c_idx]: $(setstorage[c_idx])")
                                #println("tmp_mlt: $(tmp_mlt)")
                                @inbounds pnew, xhold, xprd = mul_rev(parent_value, setstorage[c_idx], tmp_mlt)
                                #println(" output ")
                                #println(" pnew = $pnew")
                                #println(" xhold = $xhold")
                                #println(" xprd = $xprd")
                            end
                            if isempty(pnew) || isempty(xhold) || isempty(xprd)
                                continue_flag = false
                                break
                            end
                            if isnan(pnew)
                                pnew = interval_MC(parent_value)
                            end
                            if isnan(xhold)
                                xhold = interval_MC(setstorage[c_idx])
                            end
                            setstorage[k] = set_value_post(x_values,  pnew, current_node, subgrad_tighten)
                            setstorage[c_idx] = set_value_post(x_values, xhold, current_node, subgrad_tighten)
                            #println("setstorage[k] = $(setstorage[k]), setstorage[c_idx] = $(setstorage[c_idx])")
                            count += 1
                        end
                    else
                        break
                    end
                end

            elseif (op == 4) # :^

                #println("^")
                child1 = first(children_idx)
                child2 = last(children_idx)
                @assert n_children == 2
                @inbounds ix1 = children_arr[child1]
                @inbounds ix2 = children_arr[child2]
                @inbounds chdset1 = numvalued[ix1]
                @inbounds chdset2 = numvalued[ix2]
                #println("parent_value = $(parent_value), setstorage[ix1] = $(setstorage[ix1]), numberstorage[ix2] = $(numberstorage[ix2])")
                if chdset1
                    pnew, xnew, ynew = power_rev(parent_value, numberstorage[ix1], setstorage[ix2])
                elseif chdset2
                    pnew, xnew, ynew = power_rev(parent_value, setstorage[ix1], numberstorage[ix2])
                else
                    pnew, xnew, ynew = power_rev(parent_value, setstorage[ix1], setstorage[ix2])
                end
                #println("pnew, xnew, ynew = $pnew, $xnew, $ynew")
                if (isempty(pnew) || isempty(xnew) || isempty(ynew))
                    continue_flag = false
                    break
                end
                if isnan(pnew)
                    pnew = interval_MC(parent_value)
                end
                setstorage[k] = pnew
                if ~chdset1
                    if isnan(xnew)
                        xnew = interval_MC(setstorage[ix1])
                    end
                    setstorage[ix1] = set_value_post(x_values, xnew, current_node, subgrad_tighten)
                end
                if ~chdset2
                    if isnan(ynew)
                        ynew = interval_MC(setstorage[ix2])
                    end
                    setstorage[ix2] = set_value_post(x_values, ynew, current_node, subgrad_tighten)
                end

                #println("setstorage[$k] = $(setstorage[k]), setstorage[ix1] = $(setstorage[ix1]), setstorage[ix2] = $(setstorage[ix2])")
            elseif (op == 5) # :/
                #println("/")

                child1 = first(children_idx)
                child2 = last(children_idx)
                @assert n_children == 2
                @inbounds ix1 = children_arr[child1]
                @inbounds ix2 = children_arr[child2]
                @inbounds chdset1 = numvalued[ix1]
                @inbounds chdset2 = numvalued[ix2]
                if chdset1
                    pnew, xnew, ynew = div_rev(parent_value, numberstorage[ix1], setstorage[ix2])
                elseif chdset2
                    pnew, xnew, ynew = div_rev(parent_value, setstorage[ix1], numberstorage[ix2])
                else
                    pnew, xnew, ynew = div_rev(parent_value, setstorage[ix1], setstorage[ix2])
                end
                if (isempty(pnew) || isempty(xnew) || isempty(ynew))
                    continue_flag = false
                    break
                end
                setstorage[parent_index] = pnew
                if isnan(pnew)
                    pnew = interval_MC(pnew)
                end
                if ~chdset1
                    if isnan(xnew)
                        xnew = interval_MC(xnew)
                    end
                    setstorage[ix1] = set_value_post(x_values, xnew, current_node, subgrad_tighten)
                end
                if ~chdset2
                    if isnan(ynew)
                        ynew = interval_MC(ynew)
                    end
                    setstorage[ix2] = set_value_post(x_values, ynew, current_node, subgrad_tighten)
                end

            elseif (op == 6) # ifelse
                continue
            end
        elseif (nod.nodetype == JuMP._Derivatives.CALLUNIVAR) # assumes that child is set-valued and thus parent is set-valued
            op = nod.index
            child_idx = children_arr[adj.colptr[k]]
            @inbounds child_value = setstorage[child_idx]
            @inbounds parent_value = setstorage[k]
            #println("child value: $(child_value)")
            #println("parent_value: $(parent_value)")
            pnew, cnew = eval_univariate_set_reverse(op, parent_value, child_value)
            if (isempty(pnew) || isempty(cnew))
                continue_flag = false
                break
            end
            if isnan(pnew)
                pnew = interval_MC(parent_value)
            end
            if isnan(cnew)
                cnew = interval_MC(child_value)
            end
            @inbounds setstorage[child_idx] = set_value_post(x_values, cnew, current_node, subgrad_tighten)
            @inbounds setstorage[k] = set_value_post(x_values, pnew, current_node, subgrad_tighten)
            #println("setstorage[$k] = $(setstorage[k]), setstorage[$(child_idx)] = $(setstorage[child_idx])")
        end
        ~continue_flag && break
    end
    return continue_flag
end
