# perform inner check and elimination for  log(a^x) = x*log(a), ...
function flatten_log_a_x!(k, adj, y, op, children_arr)
    @inbounds child_idx = children_arr[adj.colptr[k]]
    cnod = y.nd[child_idx]
    cop = cnod.index
    if (cop == 4) && (cnod.nodetype = JuMP._Derivatives.CALL)      # Is exponential?
        @inbounds children_idx = nzrange(adj, child_idx)
        idx1 = first(children_idx)
        idx2 = last(children_idx)
        @inbounds ix1 = children_arr[idx1]     # base index
        @inbounds ix2 = children_arr[idx2]     # exponential index
        if (isnum[idx1] && ~isnum[idx2])       # Is log(a^x) = x*log(a)?
            y.nd[k] = JuMP._Derivatives.NodeData(JuMP._Derivatives.CALL, 3, y.nd[k].parent)
            y.nd[child_idx] = JuMP._Derivatives.NodeData(JuMP._Derivatives.CALLUNIVAR, op, y.nd[child_idx].parent)
            y.nd[ix1] = JuMP._Derivatives.NodeData(y.nd[ix1].nodetype, y.nd[ix1].index, child_idx)
            y.nd[ix2] = JuMP._Derivatives.NodeData(y.nd[ix2].nodetype, y.nd[ix2].index, y.nd[k].parent)
            adj[child_idx,ix2] = 0
            adj[k,ix2] = 1
            # Fix remaining node in list for appropriate modifications
        end
    end
end

# perform inner check and elimination for log(xy...z) = log(x) + log(y) ... log(z)
function flatten_log_xmy!(k, adj, y, op, children_arr)
end

# perform inner check and elimination for log(x/y) = log(x)-log(y) ...
function flatten_log_xdy!(k, adj, y, op, children_arr)
end

# perform inner check and elimination for (a^x)^b = (a^b)^x ...
function flatten_exp_axb!(k, adj, y, children_arr)
end

# perform inner check and elimination for (x^a)^b = x^(ab) ...
function flatten_exp_xab!(k, adj, y, children_arr)
    @inbounds children_idx = children_arr[adj.colptr[k]]
    child_idx1 = first(children_idx)
    child_idx2 = last(children_idx)
    cn1 = y.nd[child_idx1]
    if (~y.isnum[child_idx1] && y.isnum[child_idx2])
        if (cn1.nodetype == JuMP._Derivatives.CALL) && (cn1.index == 4)
            @inbounds subchildren_idx = children_arr[adj.colptr[child_idx1]]
            subchild_idx1 = first(subchildren_idx)
            subchild_idx2 = last(subchildren_idx)
            if (~y.isnum[subchild_idx1] && y.isnum[subchild_idx2])
                # Fix graph kernel
                y.nd[child_idx1] = JuMP._Derivatives.NodeData(y.nd[subchild_idx1].nodetype, y.nd[subchild_idx1].index, y.nd[k].parent)
                y.nd[child_idx2] = JuMP._Derivatives.NodeData(JuMP._Derivatives.CALL, 3, y.nd[k].parent)
                y.nd[subchild_idx1] = JuMP._Derivatives.NodeData(y.nd[child_idx2].nodetype, y.nd[child_idx2].index, y.nd[child_idx2])
                y.nd[subchild_idx2] = JuMP._Derivatives.NodeData(y.nd[subchild_idx2].nodetype, y.nd[subchild_idx2].index, y.nd[child_idx2])
                adj[child_idx1,subchild_idx1] = 0
                adj[child_idx1,subchild_idx2] = 0
                adj[child_idx2,subchild_idx1] = 1
                adj[child_idx2,subchild_idx2] = 1
                # Fix remaining node in list for appropriate modifications
            end
        end
    end
end

# perform inner check and elimination for a^(log(x)) = x^log(a) ...
function flatten_exp_alogx!(k, adj, y, children_arr)

    @inbounds children_idx = children_arr[adj.colptr[k]]
    child_idx1 = first(children_idx)
    child_idx2 = last(children_idx)
    cnod1 = y.nd[child_idx1]
    cnod2 = y.nd[child_idx2]

    @inbounds subchild_idx = children_arr[adj.colptr[child_idx2]]
    scnod = y.nd[subchild_idx]
    if (y.isnum[child_idx1] && ~y.isnum[child_idx2])
        if (cnod2.nodetype = JuMP._Derivatives.CALLUNIVAR)
            if ((cop2 == 8) || (cop2 == 9) || (cop2 == 10) || (cop2 == 11))
                y.nd[child_idx1] = JuMP._Derivatives.NodeData(scnod.nodetype, scnod.index, scnod.parent)
                y.nd[subchild_idx] = JuMP._Derivatives.NodeData(cnod1.nodetype, cnod1.index, cnod1.parent)
                # Fix remaining node in list for appropriate modifications
            end
        end
    end
end

# perform inner check and elimination for exp(x)exp(y) = exp(x+y) ...
function flatten_mult_expxy!(k, adj, y, children_arr)
end

"""
    flatten_function!
"""
function flatten_function!(y)
    children_arr = rowvals(adj)
    for k in length(y.nd):-1:1
        @inbounds nod = y.nd[k]
        if nod.nodetype == JuMP._Derivatives.CALL
            op = nod.index
            if (op == 3)
                flatten_mult_expxy!(k, adj, y, children_arr)      # exp(x)exp(y) = exp(x+y)
            elseif (op == 4)
                flatten_exp_axb!(k, adj, y, children_arr)         # (a^x)^b = (a^b)^x
                flatten_exp_xab!(k, adj, y, children_arr)         # (x^a)^b = x^(a*b)
                flatten_exp_alogx!(k, adj, y, children_arr)       # a^log(x) = x^log(a)
            end
        elseif nod.nodetype == JuMP._Derivatives.CALLUNIVAR
            op = nod.index
            if (op == 8) || (op == 9) || (op == 10) || (op == 11)
                flatten_log_a_x!(k, adj, y, op, children_arr)      # log(a^x) = xlog(a)
                flatten_log_xdy!(k, adj, y, op, children_arr)      # log(x/y) = log(x)-log(y)
                flatten_log_xmy!(k, adj, y, op, children_arr)      # log(xy) = log(x)+log(y)
            end
        end
    end
end

"""
    dag_flattening!
"""
function dag_flattening!(d::Evaluator)

    # flatten objective
    flatten_function!(d.objective)
    SparseArrays.dropzero!(d.objective.adj)

    # flatten constraints
    for i in length(d.constraints)
        flatten_function!(d.constraints[i])
        SparseArrays.dropzero!(d.constraints[i].adj)
    end

    # flatten subexpression
    for i in length(d.subexpressions)
        flatten_function!(d.subexpressions[i])
        SparseArrays.dropzero!(d.subexpressions[i].adj)
    end
end
