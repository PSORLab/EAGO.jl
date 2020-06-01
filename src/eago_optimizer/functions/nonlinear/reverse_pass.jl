
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
