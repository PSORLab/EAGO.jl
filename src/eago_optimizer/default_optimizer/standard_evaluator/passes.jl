"""
    set_value_construct

Construct the set-valued relaxation of variable `i` of `N` in `node` at `x_values`.
"""
function set_value_construct(i::Int, N::Int, x_values::Vector{Float64}, node::NodeBB)
    @inbounds xval = x_values[i]
    @inbounds intv = IntervalType(node.lower_variable_bounds[i], node.upper_variable_bounds[i])
    seed = seed_gradient(Float64, i, N)
    @inbounds MC{N}(xval, xval, intv, seed, seed, false)
end

"""
    set_value_post

Post process set_value operator. By default, performs the affine interval cut on
a MC structure.
"""
function set_value_post(x_values::Vector{Float64}, val::MC{N}, node::NodeBB, flag::Bool) where N
    if flag
        lower = val.cv
        upper = val.cc
        for i in 1:N
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
        return MC{N}(val.cv, val.cc, Interval(lower, upper), val.cv_grad, val.cc_grad, val.cnst)
    else
        return val
    end
end

id_to_operator = Dict(value => key for (key, value) in JuMP.univariate_operator_to_id)
function forward_eval(setstorage::Vector{MC{N}}, numberstorage::Vector{Float64}, numvalued::Vector{Bool},
                      nd::Vector{JuMP.NodeData}, adj::SparseMatrixCSC{Bool,Int}, const_values::Vector{Float64},
                      parameter_values::Vector{Float64}, current_node::NodeBB, x_values::Vector{Float64},
                      subexpr_values_flt::Vector{Float64}, subexpr_values_set::Vector{MC{N}},
                      subexpression_isnum::Vector{Bool}, user_input_buffer::Vector{MC{N}}, subgrad_tighten::Bool,
                      tpdict::Dict{Int,Tuple{Int,Int,Int,Int}}, tp1storage::Vector{Float64},
                      tp2storage::Vector{Float64}, tp3storage::Vector{Float64}, tp4storage::Vector{Float64},
                      first_eval_flag::Bool;
                      user_operators::JuMP._Derivatives.UserOperatorRegistry = JuMP._Derivatives.UserOperatorRegistry()) where N

    @assert length(numberstorage) >= length(nd)
    @assert length(setstorage) >= length(nd)
    @assert length(numvalued) >= length(nd)

    set_value_sto = zero(MC{N})

    children_arr = rowvals(adj)
    n = length(x_values)

    for k in length(nd):-1:1
        @inbounds nod = nd[k]
        op = nod.index
        if nod.nodetype == JuMP._Derivatives.VARIABLE                                   # DONE
            set_value_sto = set_value_construct(nod.index, N, x_values, current_node)
            if first_eval_flag
                setstorage[k] = set_value_sto
            else
                setstorage[k] = set_value_sto ∩ setstorage[k]
            end
            numvalued[k] = false
        elseif nod.nodetype == JuMP._Derivatives.VALUE                                  # DONE
            @inbounds numberstorage[k] = const_values[nod.index]
            numvalued[k] = true
        elseif nod.nodetype == JuMP._Derivatives.SUBEXPRESSION                          # DONE
            @inbounds isnum = subexpression_isnum[nod.index]
            if isnum
                @inbounds numberstorage[k] = subexpr_values_flt[nod.index]
            else
                @inbounds setstorage[k] = subexpr_values_set[nod.index]
            end
            numvalued[k] = isnum
        elseif nod.nodetype == JuMP._Derivatives.PARAMETER                              # DONE
            @inbounds numberstorage[k] = parameter_values[nod.index]
            numvalued[k] = true
        elseif nod.nodetype == JuMP._Derivatives.CALL
            op = nod.index
            @inbounds children_idx = nzrange(adj,k)
            n_children = length(children_idx)
            if op == 1 # :+
                tmp_sum = 0.0
                isnum = true
                chdset = true
                for c_idx in children_idx
                    @inbounds ix = children_arr[c_idx]
                    @inbounds chdset = numvalued[ix]
                    if (chdset)
                        @inbounds tmp_sum += numberstorage[ix]
                    else
                        if first_eval_flag
                            @inbounds tmp_sum += setstorage[ix]
                        else
                            @inbounds tmp_sum = McCormick.plus_kernel(tmp_sum, setstorage[ix], setstorage[k].Intv)
                        end
                    end
                    @inbounds isnum &= chdset
                end
                numvalued[k] = isnum
                if (isnum)
                    numberstorage[k] = tmp_sum
                else
                    setstorage[k] = tmp_sum
                end                                                        #DONE
            elseif op == 2 # :-                                                     # DONE
                child1 = first(children_idx)
                @assert n_children == 2
                @inbounds ix1 = children_arr[child1]
                @inbounds ix2 = children_arr[child1+1]
                @inbounds chdset1 = numvalued[ix1]
                @inbounds chdset2 = numvalued[ix2]
                @inbounds isnum = chdset1
                @inbounds isnum &= chdset2
                if chdset1
                    @inbounds tmp_sub = numberstorage[ix1]
                else
                    @inbounds tmp_sub = setstorage[ix1]
                end
                if chdset2
                    @inbounds tmp_sub -= numberstorage[ix2]
                else
                    if first_eval_flag
                        @inbounds tmp_sub -= setstorage[ix2]
                    else
                        @inbounds tmp_sub = McCormick.minus_kernel(tmp_sub, setstorage[ix2], setstorage[k].Intv)
                    end
                end
                numvalued[k] = isnum
                if (isnum)
                    numberstorage[k] = tmp_sub
                else
                    setstorage[k] = set_value_post(x_values, tmp_sub, current_node, subgrad_tighten)
                end
            elseif op == 3 # :*
                tmp_prod = 1.0
                isnum = true
                chdset = true
                for c_idx in children_idx
                    @inbounds chdset = numvalued[children_arr[c_idx]]
                    isnum &= chdset
                    if (chdset)
                        @inbounds tmp_prod *= numberstorage[children_arr[c_idx]]
                    else
                        if first_eval_flag
                            @inbounds tmp_prod = tmp_prod*setstorage[children_arr[c_idx]]
                        else
                            @inbounds tmp_prod = McCormick.mult_kernel(tmp_prod, setstorage[children_arr[c_idx]], setstorage[k].Intv)
                        end
                    end
                end
                if (isnum)
                    numberstorage[k] = tmp_prod
                else
                    setstorage[k] = set_value_post(x_values, tmp_prod, current_node, subgrad_tighten)
                end
                numvalued[k] = isnum
            elseif op == 4 # :^                                                      # DONE
                @assert n_children == 2
                idx1 = first(children_idx)
                idx2 = last(children_idx)
                @inbounds ix1 = children_arr[idx1]
                @inbounds ix2 = children_arr[idx2]
                @inbounds chdset1 = numvalued[ix1]
                @inbounds chdset2 = numvalued[ix2]
                if chdset1
                    @inbounds base = numberstorage[ix1]
                else
                    @inbounds base = setstorage[ix1]
                end
                if chdset2
                    @inbounds exponent = numberstorage[ix2]
                else
                    @inbounds exponent = setstorage[ix2]
                end
                if exponent == 1
                    if chdset1
                        @inbounds numberstorage[k] = base
                    else
                        @inbounds setstorage[k] = base
                    end
                else
                    if (chdset1 && chdset2)
                        numberstorage[k] = base^exponent
                    else
                        if first_eval_flag
                            tmp_pow = pow(base, exponent)
                        else
                            @inbounds tmp_pow = ^(base, exponent, setstorage[k].Intv)
                        end
                        setstorage[k] = set_value_post(x_values, tmp_pow, current_node, subgrad_tighten)
                    end
                end
                numvalued[k] = (chdset1 && chdset2)
            elseif op == 5 # :/                                                      # DONE
                @assert n_children == 2
                idx1 = first(children_idx)
                idx2 = last(children_idx)
                @inbounds ix1 = children_arr[idx1]
                @inbounds ix2 = children_arr[idx2]
                @inbounds chdset1 = numvalued[ix1]
                @inbounds chdset2 = numvalued[ix2]
                if chdset1
                    @inbounds numerator = numberstorage[ix1]
                else
                    @inbounds numerator = setstorage[ix1]
                end
                if chdset2
                    @inbounds denominator = numberstorage[ix2]
                else
                    @inbounds denominator = setstorage[ix2]
                end
                if chdset1 && chdset2
                    numberstorage[k] = numerator/denominator
                else
                    if first_eval_flag
                        tmp_div = numerator/denominator
                    else
                        @inbounds tmp_div = div_kernel(numerator, denominator, setstorage[k].Intv)
                    end
                    setstorage[k] = set_value_post(x_values, tmp_div, current_node, subgrad_tighten)
                end
                numvalued[k] = chdset1 && chdset2
            elseif op == 6 # ifelse
                @assert n_children == 3
                idx1 = first(children_idx)
                @inbounds chdset1 = numvalued[idx1]
                if chdset1
                    @inbounds condition = setstorage[children_arr[idx1]]
                else
                    @inbounds condition = numberstorage[children_arr[idx1]]
                end
                @inbounds chdset2 = numvalued[children_arr[idx1+1]]
                @inbounds chdset3 = numvalued[children_arr[idx1+2]]
                if chdset2
                    @inbounds lhs = setstorage[children_arr[idx1+1]]
                else
                    @inbounds lhs = numberstorage[children_arr[idx1+1]]
                end
                if chdset3
                    @inbounds rhs = setstorage[children_arr[idx1+2]]
                else
                    @inbounds rhs = numberstorage[children_arr[idx1+2]]
                end
                error("IF ELSE TO BE IMPLEMENTED SHORTLY")
                #storage[k] = set_value_post(x_values, ifelse(condition == 1, lhs, rhs), current_node)                                                 # DONE
            elseif op >= JuMP._Derivatives.USER_OPERATOR_ID_START
                op_sym = id_to_operator[op]
                evaluator = user_operators.multivariate_operator_evaluator[op - JuMP._Derivatives.USER_OPERATOR_ID_START+1]
                f_input = view(user_input_buffer, 1:n_children)
                r = 1
                isnum = true
                for c_idx in children_idx
                    @inbounds ix = children_arr[c_idx]
                    @inbounds chdset = numvalued[ix]
                    isnum &= chdset
                    if chdset
                        @inbounds f_input[r] = numberstorage[ix]
                    else
                        @inbounds f_input[r] = setstorage[ix]
                    end
                    r += 1
                end
                fval = MOI.eval_objective(evaluator, f_input)
                if isnum
                    numberstorage[k] = fval
                else
                    setstorage[k] = set_value_post(x_values, fval, current_node, subgrad_tighten)
                end
                numvalued[k] = isnum                  # DONE
            else
                error("Unsupported operation $(operators[op])")
            end
        elseif nod.nodetype == JuMP._Derivatives.CALLUNIVAR                         # DONE
            @inbounds child_idx = children_arr[adj.colptr[k]]
            @inbounds chdset = numvalued[child_idx]
            if chdset
                @inbounds child_val = numberstorage[child_idx]
            else
                @inbounds child_val = setstorage[child_idx]
            end
            if op >= JuMP._Derivatives.USER_UNIVAR_OPERATOR_ID_START
                userop = op - JuMP._Derivatives.USER_UNIVAR_OPERATOR_ID_START + 1
                @inbounds f = user_operators.univariate_operator_f[userop]
                fval = f(child_val)
                if chdset
                    @inbounds numberstorage[k] = fval
                else
                    @inbounds setstorage[k] = set_value_post(x_values, fval, current_node, subgrad_tighten)
                end
            elseif single_tp(op) && chdset
                @inbounds tpdict_tuple = tpdict[k]
                tindx1 = tpdict_tuple[1]
                tindx2 = tpdict_tuple[2]
                @inbounds tp1 = tp1storage[tindx1]
                @inbounds tp2 = tp2storage[tindx2]
                fval, tp1, tp2 = single_tp_set(op, child_val, setstorage[k], tp1, tp2, first_eval_flag)
                @inbounds tp1storage[tindx] = tp1
                @inbounds tp1storage[tindx] = tp2
                @inbounds setstorage[k] = set_value_post(x_values, fval, current_node, subgrad_tighten)
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
                fval, tp1, tp2, tp3, tp4 = double_tp_set(op, child_val, setstorage[k], tp1, tp2, tp3, tp4, first_eval_flag)
                @inbounds tp1storage[tindx1] = tp1
                @inbounds tp2storage[tindx2] = tp2
                @inbounds tp3storage[tindx1] = tp3
                @inbounds tp4storage[tindx2] = tp4
                @inbounds setstorage[k] = set_value_post(x_values, fval, current_node, subgrad_tighten)
            else
                fval = eval_univariate_set(op, child_val)
                if chdset
                    @inbounds numberstorage[k] = fval
                else
                    @inbounds setstorage[k] = set_value_post(x_values, fval, current_node, subgrad_tighten)
                end
            end
            numvalued[k] = chdset
        elseif nod.nodetype == JuMP._Derivatives.COMPARISON                         # DONE
            op = nod.index
            @inbounds children_idx = nzrange(adj,k)
            n_children = length(children_idx)
            result = true
            for r in 1:n_children-1
                @inbounds ix1 = children_arr[children_idx[r]]
                @inbounds ix2 = children_arr[children_idx[r+1]]
                @inbounds isnum1 = numvalued[ix1]
                @inbounds isnum2 = numvalued[ix2]
                if isnum1
                    @inbounds cval_lhs = numberstorage[ix1]
                else
                    @inbounds cval_lhs = setstorage[ix1]
                end
                if isnum2
                    @inbounds cval_rhs = numberstorage[ix2]
                else
                    @inbounds cval_rhs = setstorage[ix2]
                end
                if op == 1
                    result &= cval_lhs <= cval_rhs
                elseif op == 2
                    result &= cval_lhs == cval_rhs
                elseif op == 3
                    result &= cval_lhs >= cval_rhs
                elseif op == 4
                    result &= cval_lhs < cval_rhs
                elseif op == 5
                    result &= cval_lhs > cval_rhs
                end
            end
            numberstorage[k] = result
        elseif nod.nodetype == JuMP._Derivatives.LOGIC                              # DONE
            op = nod.index
            @inbounds children_idx = nzrange(adj,k)
            ix1 = children_arr[first(children_idx)]
            ix2 = children_arr[last(children_idx)]
            @inbounds isnum1 = numvalued[ix1]
            @inbounds isnum2 = numvalued[ix2]
            if isnum1
                cval_lhs = (numberstorage[ix1] == 1)
            else
                cval_lhs = (setstorage[ix1] == 1)
            end
            if isnum2
                cval_rhs = (numberstorage[ix2] == 1)
            else
                cval_rhs = (setstorage[ix2] == 1)
            end
            if op == 1
                numberstorage[k] = cval_lhs && cval_rhs
            elseif op == 2
                numberstorage[k] = cval_lhs || cval_rhs
            end
        else                                                                        # DONE
            error("Unrecognized node type $(nod.nodetype).")
        end
    end
    if numvalued[1]
        return numberstorage[1]
    else
        return setstorage[1]
    end
end

function forward_eval_obj(d::Evaluator, x::Vector{Float64})
    subexpr_values_flt = d.subexpression_values_flt
    subexpr_values_set = d.subexpression_values_set
    user_operators = d.m.nlp_data.user_operators::JuMP._Derivatives.UserOperatorRegistry
    user_input_buffer = d.jac_storage
    subgrad_tighten = d.subgrad_tighten

    for (ind, k) in enumerate(reverse(d.subexpression_order))
        ex = d.subexpressions[k]
        temp = forward_eval(ex.setstorage, ex.numberstorage, ex.numvalued,
                                         ex.nd, ex.adj, ex.const_values,
                                         d.parameter_values, d.current_node,
                                         x, subexpr_values_flt, subexpr_values_set, d.subexpression_isnum,
                                         user_input_buffer, subgrad_tighten, #ex.tpdict,
                                         user_operators=user_operators)
        d.subexpression_isnum[ind] = ex.numvalued[1]
        if d.subexpression_isnum[ind]
            d.subexpression_values_flt[k] = temp
        else
            d.subexpression_values_set[k] = temp
        end
    end

    if d.has_nlobj
        ex = d.objective
        forward_eval(ex.setstorage, ex.numberstorage, ex.numvalued,
                     ex.nd, ex.adj, ex.const_values,
                     d.parameter_values, d.current_node,
                     x, subexpr_values_flt, subexpr_values_set, d.subexpression_isnum,
                     user_input_buffer, subgrad_tighten, #ex.tpdict,
                     user_operators=user_operators)
    end
end

function forward_eval_all(d::Evaluator, x::Vector{Float64})

    subexpr_values_flt = d.subexpression_values_flt
    subexpr_values_set = d.subexpression_values_set
    user_operators = d.m.nlp_data.user_operators::JuMP._Derivatives.UserOperatorRegistry
    user_input_buffer = d.jac_storage
    subgrad_tighten = d.subgrad_tighten
    first_eval_flag = d.first_eval_flag

    for (ind, k) in enumerate(reverse(d.subexpression_order))
        ex = d.subexpressions[k]
        temp = forward_eval(ex.setstorage, ex.numberstorage, ex.numvalued,
                            ex.nd, ex.adj, ex.const_values,
                            d.parameter_values, d.current_node,
                            x, subexpr_values_flt, subexpr_values_set, d.subexpression_isnum,
                            user_input_buffer, subgrad_tighten, ex.tpdict,
                            ex.tp1storage, ex.tp2storage, ex.tp3storage, ex.tp4storage,
                            first_eval_flag, user_operators=user_operators)
        d.subexpression_isnum[ind] = ex.numvalued[1]
        if d.subexpression_isnum[ind]
            d.subexpression_values_flt[k] = temp
        else
            d.subexpression_values_set[k] = temp
        end
    end

    if d.has_nlobj
        ex = d.objective
        forward_eval(ex.setstorage, ex.numberstorage, ex.numvalued,
                     ex.nd, ex.adj, ex.const_values,
                     d.parameter_values, d.current_node,
                     x, subexpr_values_flt, subexpr_values_set,
                     d.subexpression_isnum, user_input_buffer, subgrad_tighten, ex.tpdict,
                     ex.tp1storage, ex.tp2storage, ex.tp3storage, ex.tp4storage,
                     first_eval_flag, user_operators=user_operators)
    end

    for (ind,ex) in enumerate(d.constraints)
        forward_eval(ex.setstorage, ex.numberstorage, ex.numvalued,
                     ex.nd, ex.adj, ex.const_values,
                     d.parameter_values, d.current_node,
                     x, subexpr_values_flt, subexpr_values_set,
                     d.subexpression_isnum, user_input_buffer, subgrad_tighten, ex.tpdict,
                     ex.tp1storage, ex.tp2storage, ex.tp3storage, ex.tp4storage,
                     first_eval_flag, user_operators=user_operators)
    end
end

# maximum number to perform reverse operation on associative term by summing and evaluating pairs
# remaining terms not reversed
const MAX_ASSOCIATIVE_REVERSE = 4
function reverse_eval(setstorage::Vector{T}, numberstorage, numvalued, subexpression_isnum,
                      subexpr_values_set, nd::Vector{JuMP.NodeData}, adj, x_values, current_node::NodeBB, subgrad_tighten::Bool) where T

    @assert length(setstorage) >= length(nd)
    @assert length(numberstorage) >= length(nd)
    @assert length(numvalued) >= length(nd)

    children_arr = rowvals(adj)
    N = length(x_values)

    tmp_hold = zero(T)
    continue_flag = true

    for k in 1:length(nd)
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
            @inbounds isnum = subexpression_isnum[nod.index]
            if ~isnum
                @inbounds subexpr_values_set[nod.index] = setstorage[k]
            end          # DONE
        elseif numvalued[k]
            continue                                             # DONE
        elseif (nod.nodetype == JuMP._Derivatives.CALL)
            op = nod.index
            parent_index = nod.parent
            @inbounds children_idx = nzrange(adj,k)
            @inbounds parent_value = setstorage[k]
            n_children = length(children_idx)
            if (op >= JuMP._Derivatives.USER_OPERATOR_ID_START)
                continue
            elseif (op == 1) # :+
                lenx = length(children_idx)
                if false #lenx == 2
                    #=
                    @inbounds child1 = first(children_idx)
                    @inbounds child2 = last(children_idx)
                    println("child1: $child1")
                    println("child2: $child2")
                    @inbounds ix1 = children_arr[child1]
                    @inbounds ix2 = children_arr[child2]
                    println("ix1: $ix1")
                    println("ix2: $ix2")
                    @inbounds chdset1 = numvalued[ix1]
                    @inbounds chdset2 = numvalued[ix2]
                    @inbounds nsto1 = numberstorage[ix1]
                    @inbounds nsto2 = numberstorage[ix2]
                    @inbounds ssto1 = setstorage[ix1]
                    @inbounds ssto2 = setstorage[ix2]
                    if (chdset1 && ~chdset2)
                        pnew, xhold, xsum = plus_rev(parent_value, nsto1, ssto2)
                        continue_flag &= ~isempty(xsum)
                        if continue_flag
                            @inbounds setstorage[ix2] = xsum
                        end
                    elseif (~chdset1 && chdset2)
                        pnew, xhold, xsum = plus_rev(parent_value, ssto1, nsto2)
                        continue_flag &= ~isempty(xhold)
                        if continue_flag
                            @inbounds setstorage[ix1] = xhold
                        end
                    elseif (~chdset1 && ~chdset2)
                        pnew, xhold, xsum = plus_rev(parent_value, ssto1, ssto2)
                        continue_flag &= ~isempty(xhold)
                        continue_flag &= ~isempty(xsum)
                        if continue_flag
                            @inbounds setstorage[ix1] = xhold
                            @inbounds setstorage[ix2] = xsum
                        end
                    end
                    continue_flag &= ~isempty(pnew)
                    if continue_flag
                        @inbounds setstorage[k] = pnew
                    else
                        println("infeasible +")
                        break
                    end
                    =#
                else
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
                    setstorage[k] = pnew
                    if ~chdset1
                        setstorage[ix1] = set_value_post(x_values, xnew, current_node, subgrad_tighten)
                    end
                    if  ~chdset2
                        setstorage[ix2] = set_value_post(x_values, ynew, current_node, subgrad_tighten)
                    end
                end
            elseif (op == 3) # :*
                tmp_sum = 1.0
                chdset = true
                count = 0
                child_arr_indx = children_arr[children_idx]
                for c_idx in child_arr_indx
                    if (count < MAX_ASSOCIATIVE_REVERSE)
                        if ~numvalued[c_idx]
                            tmp_sum = 1.0
                            for cin_idx in child_arr_indx
                                if (cin_idx != c_idx)
                                    @inbounds chdset = numvalued[cin_idx]
                                    if (chdset)
                                        @inbounds tmp_sum *= numberstorage[cin_idx]
                                    else
                                        @inbounds tmp_sum *= setstorage[cin_idx]
                                    end
                                end
                            end
                            @inbounds chdset = numvalued[c_idx]
                            if (chdset)
                                @inbounds pnew, xhold, xsum = mul_rev(parent_value, numberstorage[c_idx], tmp_sum)
                            else
                                @inbounds pnew, xhold, xsum = mul_rev(parent_value, setstorage[c_idx], tmp_sum)
                            end
                            if (isempty(pnew) || isempty(xhold) || isempty(xsum))
                                continue_flag = false
                                break
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
                setstorage[k] = pnew
                if ~chdset1
                    setstorage[ix1] = set_value_post(x_values, xnew, current_node, subgrad_tighten)
                end
                if ~chdset2
                    setstorage[ix2] = set_value_post(x_values, ynew, current_node, subgrad_tighten)
                end
            elseif (op == 5) # :/
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
                if ~chdset1
                    setstorage[ix1] = set_value_post(x_values, xnew, current_node, subgrad_tighten)
                end
                if ~chdset2
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
            pnew, cnew = eval_univariate_set_reverse(op, parent_value, child_value)
            if (isempty(pnew) || isempty(cnew))
                continue_flag = false
                break
            end
            @inbounds setstorage[k] = set_value_post(x_values, cnew, current_node, subgrad_tighten)
            @inbounds setstorage[nod.parent] = set_value_post(x_values, pnew, current_node, subgrad_tighten)
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
        ex.setstorage[1] = ex.setstorage[1] ∩ IntervalType(-Inf, d.objective_ubd)
        feas &= reverse_eval(ex.setstorage, ex.numberstorage, ex.numvalued, subexpr_isnum,
                              subexpr_values_set, ex.nd, ex.adj, x, d.current_node, subgrad_tighten)
    end
    for i in 1:length(d.constraints)
        # Cut constraints on constraint bounds & reverse
        if feas
            ex = d.constraints[i]
            println("ex.setstorage: $(ex.setstorage)")
            println("d.constraints_lbd: $(d.constraints_lbd)")
            println("d.constraints_ubd: $(d.constraints_ubd)")
            ex.setstorage[1] = ex.setstorage[1] ∩ IntervalType(d.constraints_lbd[i], d.constraints_ubd[i])
            feas &= reverse_eval(ex.setstorage, ex.numberstorage, ex.numvalued, subexpr_isnum,
                                  subexpr_values_set, ex.nd, ex.adj, x, d.current_node, subgrad_tighten)
        else
            break
        end
    end
    for k in 1:length(d.subexpression_order)
        if feas
            ex = d.subexpressions[d.subexpression_order[k]]
            ex.setstorage[1] = subexpression_values_set[d.subexpression_order[k]]
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
    forward_reverse_pass(d::Evaluator, x::Vector{Float64})

Performs a `d.fw_repeats` forward passes of the set-value evaluator each followed
by a reverse pass if `d.has_reverse` as long as the node between passes differs
by more that `d.fw_atol` at each iteration.
"""
function forward_reverse_pass(d::Evaluator, x::Vector{Float64})
    flag = true
    #if ~same_box(d.current_node, d.last_node, 0.0)
    #   d.last_node = d.current_node
        if (d.last_x != x)
            if d.has_reverse
                for i in 1:d.cp_reptitions
                    d.first_eval_flag = (i == 1)
                    if flag
                        forward_eval_all(d,x)
                        flag = reverse_eval_all(d,x)
                        ~flag && break
                        same_box(d.current_node, get_node(d), d.cp_tolerance) && break
                    end
                end
            else
                d.first_eval_flag = true
                forward_eval_all(d,x)
            end
        end
     #end
     return flag
end
