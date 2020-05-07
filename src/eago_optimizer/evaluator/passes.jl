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
function set_value_post(x_values::Vector{Float64}, val::MC{N,T}, node::NodeBB, flag::Bool) where {N, T<:RelaxTag}
    if flag
        lower = val.cv
        upper = val.cc
        for i = 1:N
            @inbounds cv_val = val.cv_grad[i]
            @inbounds cc_val = val.cc_grad[i]
            if (cv_val > 0.0)
                @inbounds lower += cv_val*(node.lower_variable_bounds[i] - x_values[i])
            else
                @inbounds lower += cv_val*(node.upper_variable_bounds[i] - x_values[i])
            end
            if (cc_val > 0.0)
                @inbounds upper += cc_val*(node.upper_variable_bounds[i] - x_values[i])
            else
                @inbounds upper += cc_val*(node.lower_variable_bounds[i] - x_values[i])
            end
        end
        lower = max(lower, val.Intv.lo)
        upper = min(upper, val.Intv.hi)
        return MC{N,T}(val.cv, val.cc, Interval{Float64}(lower, upper), val.cv_grad, val.cc_grad, val.cnst)
    else
        return val
    end
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
            if (~chdset1 & ~chdset2)
                tmp_mc = McCormick.plus_kernel(setstorage[ix1], setstorage[ix2], setstorage[k])
            elseif chdset1 & ~chdset2
                tmp_mc = McCormick.plus_kernel(numberstorage[ix1], setstorage[ix2], setstorage[k])
            elseif chdset2 & ~chdset1
                tmp_mc = McCormick.plus_kernel(setstorage[ix1], numberstorage[ix2], setstorage[k])
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
        if first_eval_flag
            @inbounds setstorage[k] = tmp_num + tmp_mc
        else
            @inbounds setstorage[k] = (tmp_num + tmp_mc) ∩ setstorage[k]
        end
    end

    return
end

function forward_minus!(k::Int64, x_values::Vector{Float64}, ix1::Int64, ix2::Int64, numvalued::Vector{Bool},
                        numberstorage::Vector{Float64}, setstorage::Vector{MC{N,T}},
                        current_node::NodeBB, subgrad_tighten::Bool, first_eval_flag::Bool) where {N, T<:RelaxTag}
    tmp_num_1 = 0.0
    tmp_mc_1 = zero(MC{N,T})
    @inbounds chdset1 = numvalued[ix1]
    @inbounds chdset2 = numvalued[ix2]
    isnum = chdset1
    isnum &= chdset2
    if first_eval_flag
        if isnum
            @inbounds tmp_num_1 = numberstorage[ix1] - numberstorage[ix2]
            @inbounds numberstorage[k] = tmp_num_1
        elseif ~chdset1 && chdset2
            @inbounds tmp_mc_1 = setstorage[ix1] - numberstorage[ix2]
        elseif chdset1 && ~chdset2
            @inbounds tmp_mc_1 = numberstorage[ix1] - setstorage[ix2]
        else
            @inbounds tmp_mc_1 = setstorage[ix1] - setstorage[ix2]
        end
    else
        if ~chdset1 && chdset2
            @inbounds tmp_mc_1 = McCormick.minus_kernel(setstorage[ix1], numberstorage[ix2], setstorage[k].Intv)
        elseif chdset1 && ~chdset2
            @inbounds tmp_mc_1 = McCormick.minus_kernel(numberstorage[ix1], setstorage[ix2], setstorage[k].Intv)
        else
            @inbounds tmp_mc_1 = McCormick.minus_kernel(setstorage[ix1], setstorage[ix2], setstorage[k].Intv)
        end
    end
    @inbounds numvalued[k] = isnum
    if ~isnum
        if first_eval_flag
            @inbounds setstorage[k] = set_value_post(x_values, tmp_mc_1, current_node, subgrad_tighten)
        else
            @inbounds setstorage[k] = set_value_post(x_values, tmp_mc_1 ∩ setstorage[k], current_node, subgrad_tighten)
        end
    end
    return
end

function forward_multiply!(k::Int64, x_values::Vector{Float64}, children_idx::UnitRange{Int64}, children_arr::Vector{Int64},
                           numvalued::Vector{Bool}, numberstorage::Vector{Float64}, setstorage::Vector{MC{N,T}},
                           current_node::NodeBB, subgrad_tighten::Bool, first_eval_flag::Bool) where {N, T<:RelaxTag}
    tmp_num_1 = 1.0
    tmp_mc_1 = one(MC{N,T})
    isnum = true
    chdset = true
    for c_idx in children_idx
        @inbounds cix = children_arr[c_idx]
        @inbounds chdset = numvalued[cix]
        isnum &= chdset
        if chdset
            @inbounds tmp_num_1 *= numberstorage[cix]
        else
            if first_eval_flag
                @inbounds tmp_mc_1 = tmp_mc_1*setstorage[cix]
            else
                @inbounds tmp_mc_1 = McCormick.mult_kernel(tmp_mc_1, setstorage[cix], setstorage[k].Intv)
            end
        end
    end
    if isnum
        numberstorage[k] = tmp_num_1
    else
        tmp_mc_1 *= tmp_num_1
        if first_eval_flag
            setstorage[k] = set_value_post(x_values, tmp_mc_1, current_node, subgrad_tighten)
        else
            setstorage[k] = set_value_post(x_values, tmp_mc_1 ∩ setstorage[k], current_node, subgrad_tighten)
        end
    end
    numvalued[k] = isnum
    return
end

function forward_power!(k::Int64, x_values::Vector{Float64}, children_idx::UnitRange{Int64},
                        children_arr::Vector{Int64}, numvalued::Vector{Bool}, numberstorage::Vector{Float64},
                        setstorage::Vector{MC{N,T}}, current_node::NodeBB, subgrad_tighten::Bool,
                        first_eval_flag::Bool,
                        ctx::GuardCtx) where {N, T<:RelaxTag}
    tmp_num_1 = 0.0
    tmp_num_2 = 0.0
    tmp_mc_1 = zero(MC{N,T})
    tmp_mc_2= zero(MC{N,T})
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
    if ((tmp_num_2 == 1.0) & chdset2)
        if chdset1
            @inbounds numberstorage[k] = tmp_num_1
        else
            @inbounds setstorage[k] = tmp_mc_1
        end
    else
        if (chdset1 & chdset2)
            @inbounds numberstorage[k] = tmp_num_1^tmp_num_2
        else
            if first_eval_flag
                if (~chdset1 & chdset2)
                    tmp_mc_1 = pow(tmp_mc_1, tmp_num_2)
                elseif (chdset1 & ~chdset2)
                    tmp_mc_1 = Cassette.overdub(ctx, pow, tmp_num_1, tmp_mc_2)
                elseif (~chdset1 & ~chdset2)
                    tmp_mc_1 = Cassette.overdub(ctx, pow, tmp_mc_1, tmp_mc_2)
                end
            else
                if (~chdset1 & chdset2)
                    @inbounds tmp_mc_1 = ^(tmp_mc_1, tmp_num_2, setstorage[k].Intv)
                elseif (chdset1 & ~chdset2)
                    @inbounds tmp_mc_1 = Cassette.overdub(ctx, ^, tmp_num_1, tmp_mc_2, setstorage[k].Intv)
                elseif (~chdset1 & ~chdset2)
                    @inbounds tmp_mc_1 = Cassette.overdub(ctx, ^, tmp_mc_1, tmp_mc_2, setstorage[k].Intv)
                end
            end
            if first_eval_flag
                @inbounds setstorage[k] = set_value_post(x_values, tmp_mc_1, current_node, subgrad_tighten)
            else
                @inbounds setstorage[k] = set_value_post(x_values, tmp_mc_1 ∩ setstorage[k], current_node, subgrad_tighten)
            end
        end
    end
    @inbounds numvalued[k] = (chdset1 & chdset2)
    return
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
            if ~chdset1 & chdset2
                tmp_mc_1 = tmp_mc_1/tmp_num_2
            elseif chdset1 & ~chdset2
                tmp_mc_1 = tmp_num_1/tmp_mc_2
            elseif ~chdset1 & ~chdset2
                tmp_mc_1 = tmp_mc_1/tmp_mc_2 #TODO: Overdub this with guard context
            end
        else
            @inbounds intv_k = setstorage[k].Intv
            if ~chdset1 & chdset2
                tmp_mc_1 = McCormick.div_kernel(tmp_mc_1, tmp_num_2, intv_k)
            elseif chdset1 & ~chdset2
                tmp_mc_1 = McCormick.div_kernel(tmp_num_1, tmp_mc_2, intv_k)
            elseif ~chdset1 & ~chdset2
                tmp_mc_1 = McCormick.div_kernel(tmp_mc_1, tmp_mc_2, intv_k) #TODO: Overdub this with guard context
            end
        end
        if first_eval_flag
            @inbounds setstorage[k] = set_value_post(x_values, tmp_mc_1, current_node, subgrad_tighten)
        else
            @inbounds setstorage[k] = set_value_post(x_values, tmp_mc_1 ∩ setstorage[k], current_node, subgrad_tighten)
        end
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
                forward_plus!(k, children_idx, children_arr, numvalued, numberstorage, setstorage, first_eval_flag)
            elseif op == 2 # :-
                @assert n_children == 2
                child1 = first(children_idx)
                @inbounds ix1 = children_arr[child1]
                @inbounds ix2 = children_arr[child1+1]
                forward_minus!(k, x_values, ix1, ix2, numvalued, numberstorage,
                               setstorage, current_node, subgrad_tighten, first_eval_flag)
            elseif op == 3 # :*
                forward_multiply!(k, x_values, children_idx, children_arr, numvalued, numberstorage,
                                  setstorage, current_node, subgrad_tighten, first_eval_flag)
            elseif op == 4 # :^                                                      # DONE
                @assert n_children == 2
                forward_power!(k, x_values, children_idx, children_arr, numvalued, numberstorage,
                               setstorage, current_node, subgrad_tighten, first_eval_flag, ctx)
            elseif op == 5 # :/                                                      # DONE
                @assert n_children == 2
                forward_divide!(k, x_values, children_idx, children_arr, numvalued, numberstorage,
                               setstorage, current_node, subgrad_tighten, first_eval_flag, ctx)
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
                    fval = Cassette.overdub(ctx, MOI.eval_objective, evaluator, f_input)
                    if first_eval_flag
                        setstorage[k] = set_value_post(x_values, fval, current_node, subgrad_tighten)
                    else
                        setstorage[k] = set_value_post(x_values, fval ∩ setstorage[k], current_node, subgrad_tighten)
                    end
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
                    fval_mc = Cassette.overdub(ctx, f, child_val_mc)
                    if first_eval_flag
                        @inbounds setstorage[k] = set_value_post(x_values, fval_mc, current_node, subgrad_tighten)
                    else
                        @inbounds setstorage[k] = set_value_post(x_values, fval_mc ∩ setstorage[k], current_node, subgrad_tighten)
                    end
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
                if first_eval_flag
                    @inbounds setstorage[k] = set_value_post(x_values, fval_mc, current_node, subgrad_tighten)
                else
                    @inbounds setstorage[k] = set_value_post(x_values, fval_mc ∩ setstorage[k], current_node, subgrad_tighten)
                end
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
                if first_eval_flag
                    @inbounds setstorage[k] = set_value_post(x_values, fval_mc, current_node, subgrad_tighten)
                else
                    @inbounds setstorage[k] = set_value_post(x_values, fval_mc ∩ setstorage[k], current_node, subgrad_tighten)
                end
            else
                if chdset
                    fval_num = eval_univariate_set(op, child_val_num)
                    @inbounds numberstorage[k] = fval_num
                else
                    fval_mc = Cassette.overdub(ctx, eval_univariate_set, op, child_val_mc)
                    if first_eval_flag
                        @inbounds setstorage[k] = set_value_post(x_values, fval_mc, current_node, subgrad_tighten)
                    else
                        @inbounds setstorage[k] = set_value_post(x_values, fval_mc ∩ setstorage[k], current_node, subgrad_tighten)
                    end
                end
            end
            numvalued[k] = chdset
                else
            error("Unrecognized node type $(nod.nodetype).")
        end
    end
    return
end

function forward_eval_obj(d::Evaluator, x::Vector{Float64})
    subexpr_values_flt = d.subexpression_values_flt
    subexpr_values_set = d.subexpression_values_set
    user_operators = d.user_operators
    user_input_buffer = d.jac_storage
    flt_user_input_buffer = d.flt_jac_storage
    subgrad_tighten = d.subgrad_tighten
    first_eval_flag = d.first_eval_flag
    seeds = d.seeds

    seeds = d.seeds
    for (ind, k) in enumerate(reverse(d.subexpression_order))
        subex = d.subexpressions[k]
        temp = forward_eval(subex.setstorage, subex.numberstorage, subex.numvalued,
                    subex.nd, subex.adj, d.current_node,
                    x, subexpr_values_flt, subexpr_values_set, d.subexpression_isnum,
                    user_input_buffer, flt_user_input_buffer, subgrad_tighten, subex.tpdict,
                    subex.tp1storage, subex.tp2storage, subex.tp3storage,
                    subex.tp4storage, first_eval_flag, user_operators, seeds, d.ctx)
        d.subexpression_isnum[ind] = subex.numvalued[1]
        if d.subexpression_isnum[ind]
            d.subexpression_values_flt[k] = temp
        else
            d.subexpression_values_set[k] = temp
        end
    end

    if d.has_nlobj
        ex = d.objective
        forward_eval(ex.setstorage, ex.numberstorage, ex.numvalued,
                     ex.nd, ex.adj, d.current_node,
                     x, subexpr_values_flt, subexpr_values_set,
                     d.subexpression_isnum, user_input_buffer, flt_user_input_buffer, subgrad_tighten, ex.tpdict,
                     ex.tp1storage, ex.tp2storage, ex.tp3storage, ex.tp4storage,
                     first_eval_flag, user_operators, seeds, d.ctx)
    end
    return
end

"""
$(TYPEDSIGNATURES)
"""
function forward_eval_all(d::Evaluator, x::Vector{Float64})

    subexpr_values_flt = d.subexpression_values_flt
    subexpr_values_set = d.subexpression_values_set
    user_operators = d.user_operators
    user_input_buffer = d.jac_storage
    flt_user_input_buffer = d.flt_jac_storage
    subgrad_tighten = d.subgrad_tighten
    first_eval_flag = d.first_eval_flag
    seeds = d.seeds

    for (ind, k) in enumerate(reverse(d.subexpression_order))
        subex = d.subexpressions[k]
        forward_eval(subex.setstorage, subex.numberstorage, subex.numvalued,
                    subex.nd, subex.adj, d.current_node,
                    x, subexpr_values_flt, subexpr_values_set, d.subexpression_isnum,
                    user_input_buffer, flt_user_input_buffer, subgrad_tighten, subex.tpdict,
                    subex.tp1storage, subex.tp2storage, subex.tp3storage,
                    subex.tp4storage, first_eval_flag, user_operators, seeds, d.ctx)

        d.subexpression_isnum[ind] = subex.numvalued[1]
        if d.subexpression_isnum[ind]
            d.subexpression_values_flt[k] = subex.numberstorage[1]
        else
            d.subexpression_values_set[k] = subex.setstorage[1]
        end
    end

    if d.has_nlobj
        ex = d.objective
        forward_eval(ex.setstorage, ex.numberstorage, ex.numvalued,
                     ex.nd, ex.adj, d.current_node,
                     x, subexpr_values_flt, subexpr_values_set,
                     d.subexpression_isnum, user_input_buffer, flt_user_input_buffer, subgrad_tighten, ex.tpdict,
                     ex.tp1storage, ex.tp2storage, ex.tp3storage, ex.tp4storage,
                     first_eval_flag, user_operators, seeds, d.ctx)
    end

    for (ind,ex) in enumerate(d.constraints)
        forward_eval(ex.setstorage, ex.numberstorage, ex.numvalued,
                     ex.nd, ex.adj, d.current_node,
                     x, subexpr_values_flt, subexpr_values_set,
                     d.subexpression_isnum, user_input_buffer, flt_user_input_buffer, subgrad_tighten, ex.tpdict,
                     ex.tp1storage, ex.tp2storage, ex.tp3storage, ex.tp4storage,
                     first_eval_flag, user_operators, seeds, d.ctx)
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
        elseif (nod.nodetype == JuMP._Derivatives.CALL)
            op = nod.index
            parent_index = nod.parent
            @inbounds children_idx = nzrange(adj,k)
            @inbounds parent_value = setstorage[k]
            n_children = length(children_idx)
            if (op >= JuMP._Derivatives.USER_OPERATOR_ID_START)
                # SKIPS USER DEFINE OPERATORS NOT BRIDGED INTO JuMP Tree
                # Representation
                continue
            elseif (op == 1) # :+
                lenx = length(children_idx)
                count = 0
                child_arr_indx = children_arr[children_idx]
                for c_idx in child_arr_indx
                    tmp_sum = 0.0
                    @inbounds inner_chdset = numvalued[c_idx]
                    if ~inner_chdset
                        if (count < MAX_ASSOCIATIVE_REVERSE)
                            for cin_idx in child_arr_indx
                                if (cin_idx != c_idx)
                                    @inbounds nchdset = numvalued[cin_idx]
                                    if (nchdset)
                                        @inbounds tmp_sum += numberstorage[cin_idx]
                                    else
                                        @inbounds tmp_sum += setstorage[cin_idx]
                                    end
                                end
                            end
                            @inbounds tmp_hold = setstorage[c_idx]
                            pnew, xhold, xsum = plus_rev(parent_value, tmp_hold, tmp_sum)
                            if (isempty(pnew) || isempty(xhold) || isempty(xsum))
                                continue_flag = false
                                break
                            end
                            if isnan(pnew)
                                pnew = interval_MC(pnew)
                            end
                            if isnan(xhold)
                                pnew = interval_MC(xhold)
                            end
                            setstorage[k] = set_value_post(x_values, pnew, current_node, subgrad_tighten)
                            setstorage[c_idx] = set_value_post(x_values, xhold, current_node, subgrad_tighten)
                        else
                            break
                        end
                    end
                    count += 1
                end
                if ~continue_flag
                    break
                end
            elseif (op == 2) # :-
                child1 = first(children_idx)
                child2 = last(children_idx)
                @assert n_children == 2
                @inbounds ix1 = children_arr[child1]
                @inbounds ix2 = children_arr[child2]
                @inbounds chdset1 = numvalued[ix1]
                @inbounds chdset2 = numvalued[ix2]
                if chdset1 && ~chdset2
                    pnew, xnew, ynew = minus_rev(parent_value, numberstorage[ix1], setstorage[ix2])
                elseif ~chdset1 && chdset2
                    pnew, xnew, ynew = minus_rev(parent_value, setstorage[ix1], numberstorage[ix2])
                elseif chdset1 && chdset2
                    pnew, xnew, ynew = minus_rev(parent_value, setstorage[ix1], setstorage[ix2])
                end
                if ~(~chdset1 && ~chdset2)
                    flagp = isempty(pnew)
                    flagx = isempty(xnew)
                    flagy = isempty(ynew)
                    if (flagp || flagx || flagy)
                        continue_flag = false
                        break
                    end
                    if isnan(pnew)
                        pnew = interval_MC(pnew)
                    end
                    setstorage[k] = pnew
                    if ~chdset1
                        if isnan(xnew)
                            xnew = interval_MC(xnew)
                        end
                        setstorage[ix1] = set_value_post(x_values, xnew, current_node, subgrad_tighten)
                    end
                    if  ~chdset2
                        if isnan(ynew)
                            ynew = interval_MC(ynew)
                        end
                        setstorage[ix2] = set_value_post(x_values, ynew, current_node, subgrad_tighten)
                    end
                end
            elseif (op == 3) # :*
                tmp_mlt = 1.0
                chdset = true
                count = 0
                child_arr_indx = children_arr[children_idx]
                for c_idx in child_arr_indx
                    if (count < MAX_ASSOCIATIVE_REVERSE)
                        if ~numvalued[c_idx]
                            tmp_mlt = 1.0
                            for cin_idx in child_arr_indx
                                if (cin_idx != c_idx)
                                    @inbounds chdset = numvalued[cin_idx]
                                    if (chdset)
                                        @inbounds tmp_mlt *= numberstorage[cin_idx]
                                    else
                                        @inbounds tmp_mlt *= setstorage[cin_idx]
                                    end
                                end
                            end
                            @inbounds chdset = numvalued[c_idx]
                            if (chdset)
                                @inbounds pnew, xhold, xprd = mul_rev(parent_value, numberstorage[c_idx], tmp_mlt)
                            else
                                @inbounds pnew, xhold, xprd = mul_rev(parent_value, setstorage[c_idx], tmp_mlt)
                            end
                            if (isempty(pnew) || isempty(xhold) || isempty(xprd))
                                continue_flag = false
                                break
                            end
                            if isnan(pnew)
                                pnew = interval_MC(pnew)
                            end
                            if isnan(xhold)
                                xhold = interval_MC(xhold)
                            end
                            setstorage[k] = set_value_post(x_values,  pnew, current_node, subgrad_tighten)
                            setstorage[c_idx] = set_value_post(x_values, xhold, current_node, subgrad_tighten)
                            count += 1
                        end
                    else
                        break
                    end
                end
            elseif (op == 4) # :^
                child1 = first(children_idx)
                child2 = last(children_idx)
                @assert n_children == 2
                @inbounds ix1 = children_arr[child1]
                @inbounds ix2 = children_arr[child2]
                @inbounds chdset1 = numvalued[ix1]
                @inbounds chdset2 = numvalued[ix2]
                if chdset1
                    pnew, xnew, ynew = power_rev(parent_value, numberstorage[ix1], setstorage[ix2])
                elseif chdset2
                    pnew, xnew, ynew = power_rev(parent_value, setstorage[ix1], numberstorage[ix2])
                else
                    pnew, xnew, ynew = power_rev(parent_value, setstorage[ix1], setstorage[ix2])
                end
                if (isempty(pnew) || isempty(xnew) || isempty(ynew))
                    continue_flag = false
                    break
                end
                if isnan(pnew)
                    pnew = interval_MC(pnew)
                end
                setstorage[k] = pnew
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
            elseif (op == 5) # :/
                #=
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
                =#
            elseif (op == 6) # ifelse
                continue
            end
        elseif (nod.nodetype == JuMP._Derivatives.CALLUNIVAR) # assumes that child is set-valued and thus parent is set-valued
            op = nod.index
            child_idx = children_arr[adj.colptr[k]]
            @inbounds child_value = setstorage[child_idx]
            @inbounds parent_value = setstorage[k]
            pnew, cnew = eval_univariate_set_reverse(op, parent_value, child_value)
            if (isempty(pnew) || isempty(cnew))
                continue_flag = false
                break
            end
            if isnan(pnew)
                pnew = interval_MC(pnew)
            end
            if isnan(cnew)
                cnew = interval_MC(cnew)
            end
            @inbounds setstorage[child_idx] = set_value_post(x_values, cnew, current_node, subgrad_tighten)
            @inbounds setstorage[k] = set_value_post(x_values, pnew, current_node, subgrad_tighten)
        end
        ~continue_flag && break
    end
    return continue_flag
end

# looks good
function reverse_eval_all(d::Evaluator, x::Vector{Float64})
    subexpr_values_set = d.subexpression_values_set
    subexpr_isnum = d.subexpression_isnum
    feas = true
    subgrad_tighten = d.subgrad_tighten

    if d.has_nlobj
        # Cut Objective at upper bound
        ex = d.objective
        reverse_updated_mc = ex.setstorage[1] ∩ Interval{Float64}(-Inf, d.objective_ubd)
        if isnan(reverse_updated_mc)
            reverse_updated_mc = interval_MC(reverse_updated_mc)
        end
        ex.setstorage[1] = reverse_updated_mc
        feas &= reverse_eval(ex.setstorage, ex.numberstorage, ex.numvalued, subexpr_isnum,
                             subexpr_values_set, ex.nd, ex.adj, x, d.current_node, subgrad_tighten)
    end
    for i = 1:length(d.constraints)
        # Cut constraints on constraint bounds & reverse
        if feas
            ex = d.constraints[i]
            reverse_updated_mc = ex.setstorage[1] ∩ Interval{Float64}(d.constraints_lbd[i], d.constraints_ubd[i])
            if isnan(reverse_updated_mc)
                reverse_updated_mc = interval_MC(reverse_updated_mc)
            end
            ex.setstorage[1] = reverse_updated_mc
            feas &= reverse_eval(ex.setstorage, ex.numberstorage, ex.numvalued, subexpr_isnum,
                                subexpr_values_set, ex.nd, ex.adj, x, d.current_node, subgrad_tighten)
        else
            break
        end
    end
    for k = 1:length(d.subexpression_order)
        if feas
            ex = d.subexpressions[d.subexpression_order[k]]
            ex.setstorage[1] = subexpr_values_set[d.subexpression_order[k]]
            feas &= reverse_eval(ex.setstorage, ex.numberstorage, ex.numvalued, subexpr_isnum,
                                  subexpr_values_set, ex.nd, ex.adj, x, d.current_node, subgrad_tighten)
        else
            break
        end
    end
    copyto!(d.last_x, x)

    return feas
end

"""
$(FUNCTIONNAME)

Performs a `d.fw_repeats` forward passes of the set-value evaluator each followed
by a reverse pass if `d.has_reverse` as long as the node between passes differs
by more that `d.fw_atol` at each iteration.
"""
function forward_reverse_pass(d::Evaluator, x::Vector{Float64})
    flag = true
    converged_flag = false
    if (d.last_x != x)
        if d.has_reverse
            for i in 1:d.cp_repetitions
                d.first_eval_flag = (i == 1)
                if flag
                    forward_eval_all(d,x)
                    flag = reverse_eval_all(d,x)
                    ~flag && break
                    converged_flag = same_box(d.current_node, get_node(d), d.cp_tolerance)
                    converged_flag && break
                end
            end
            if flag
                forward_eval_all(d,x)
            end
        else
            d.first_eval_flag = true
            forward_eval_all(d,x)
        end
    end
    d.last_x .= x

     return flag
end
