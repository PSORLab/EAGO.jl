
# maximum number to perform reverse operation on associative term by summing
# and evaluating pairs remaining terms not reversed
const MAX_ASSOCIATIVE_REVERSE = 4

function reverse_plus_binary!()
    return nothing
end

function reverse_plus_narity!()
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
                pnew, xhold, xsum = plus_rev(parent_value, tmp_hold, tmp_sum)
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
    return nothing
end

function reverse_multiply_binary!()
    return nothing
end

function reverse_multiply_narity!()
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
                        if chdset
                            @inbounds tmp_mlt *= numberstorage[cin_idx]
                        else
                            @inbounds tmp_mlt *= setstorage[cin_idx]
                        end
                    end
                end
                @inbounds chdset = numvalued[c_idx]
                if chdset
                    @inbounds pnew, xhold, xprd = mul_rev(parent_value, numberstorage[c_idx], tmp_mlt)
                else
                    @inbounds pnew, xhold, xprd = mul_rev(parent_value, setstorage[c_idx], tmp_mlt)
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
                count += 1
            end
        else
            break
        end
    end
    return nothing
end

function reverse_minus!()
    child1 = first(children_idx)
    child2 = last(children_idx)
    @assert n_children == 2
    @inbounds ix1 = children_arr[child1]
    @inbounds ix2 = children_arr[child2]
    @inbounds chdset1 = numvalued[ix1]
    @inbounds chdset2 = numvalued[ix2]
    if chdset1 && !chdset2
        pnew, xnew, ynew = minus_rev(parent_value, numberstorage[ix1], setstorage[ix2])
    elseif !chdset1 && chdset2
        pnew, xnew, ynew = minus_rev(parent_value, setstorage[ix1], numberstorage[ix2])
    elseif chdset1 && chdset2
        pnew, xnew, ynew = minus_rev(parent_value, setstorage[ix1], setstorage[ix2])
    end
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
            setstorage[ix1] = set_value_post(x_values, xnew, current_node, subgrad_tighten)
        end
        if  ~chdset2
            if isnan(ynew)
                ynew = interval_MC(setstorage[ix2])
            end
            setstorage[ix2] = set_value_post(x_values, ynew, current_node, subgrad_tighten)
        end
    end
    return nothing
end

function reverse_power!()
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

    return nothing
end

function reverse_divide!()
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

    return nothing
end

function reverse_univariate!()
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
        pnew = interval_MC(parent_value)
    end
    if isnan(cnew)
        cnew = interval_MC(child_value)
    end
    @inbounds setstorage[child_idx] = set_value_post(x_values, cnew, current_node, subgrad_tighten)
    @inbounds setstorage[k] = set_value_post(x_values, pnew, current_node, subgrad_tighten)

    return nothing
end

"""
$(TYPEDSIGNATURES)
"""
function reverse_pass_kernel(setstorage::Vector{T}, numberstorage, numvalued, subexpression_isnum,
                             subexpr_values_set, nd::Vector{JuMP.NodeData}, adj, x_values, current_node::NodeBB,
                             subgrad_tighten::Bool) where T

    @assert length(setstorage) >= length(nd)
    @assert length(numberstorage) >= length(nd)
    @assert length(numvalued) >= length(nd)

    children_arr = rowvals(adj)
    N = length(x_values)

    tmp_hold = zero(T)
    continue_flag = true

    for k = 1:length(nd)

        @inbounds nod = nd[k]
        ntype = nod.nodetype
        nvalued = @inbounds numvalued[k]

        if ntype == JuMP._Derivatives.VALUE      || ntype == JuMP._Derivatives.LOGIC     ||
           ntype == JuMP._Derivatives.COMPARISON || ntype == JuMP._Derivatives.PARAMETER ||
           ntype == JuMP._Derivatives.EXTRA
           continue

        elseif nod.nodetype == JuMP._Derivatives.VARIABLE
            op = nod.index
            @inbounds current_node.lower_variable_bounds[op] = setstorage[k].Intv.lo
            @inbounds current_node.upper_variable_bounds[op] = setstorage[k].Intv.hi

        elseif nod.nodetype == JuMP._Derivatives.SUBEXPRESSION
            #=
            @inbounds isnum = subexpression_isnum[nod.index]
            if ~isnum
                @inbounds subexpr_values_set[nod.index] = setstorage[k]
            end          # DONE
            =#

        elseif nvalued
            continue

        elseif nod.nodetype == JuMP._Derivatives.CALL
            op = nod.index
            parent_index = nod.parent
            @inbounds children_idx = nzrange(adj,k)
            @inbounds parent_value = setstorage[k]
            n_children = length(children_idx)

            # SKIPS USER DEFINE OPERATORS NOT BRIDGED INTO JuMP Tree Representation
            if op >= JuMP._Derivatives.USER_OPERATOR_ID_START
                continue

            # :+
            elseif op === 1
                if n_children === 2
                    reverse_plus_binary!()
                else
                    reverse_plus_narity!()
                end

            # :-
            elseif op === 2
                reverse_minus!()

            elseif op === 3 # :*
                if n_children === 2
                    reverse_multiply_binary!()
                else
                    reverse_multiply_narity!()
                end

             # :^
            elseif op === 4
                reverse_power!()

            # :/
            elseif op === 5
                reverse_divide!()

            elseif op == 6 # ifelse
                continue
            end

        # assumes that child is set-valued and thus parent is set-valued (since isnumber already checked)
        elseif nod.nodetype == JuMP._Derivatives.CALLUNIVAR
            reverse_univariate!()

        end

        !continue_flag && break
    end

    return continue_flag
end
