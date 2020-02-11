function tape_to_list(tape::Tape)
    len = length(tape.nd)
    @inbounds last_node = tape.nd[len]
    new_nds = NodeData[NodeData(last_node.nodetype, last_node.index, -1)]

    queue = Tuple{Int,Int}[(len, -1)]
    parent_dict = Dict{Int,Int}(len => 1) # starting node is 1
    node_count = 1

    while ~isempty(queue)
        println("pre queue: $(queue)")
        (active_node_num, prior_prt) = popfirst!(queue)
        @inbounds active_node = tape.nd[active_node_num]
        @inbounds active_node_child1 = active_node.children[1]
        println("CHILDREN OF ($active_node_num, $(active_node.children), $(active_node.index)) with parent = $prior_prt")
        if (active_node_child1 > 0) # has any children
            for child in active_node.children
                push!(queue, (child, active_node_num))
                @inbounds cn = tape.nd[child]
                @inbounds parent_num = parent_dict[active_node_num]
                added_node = NodeData(cn.nodetype, cn.index, parent_num)
                push!(new_nds, added_node)
                node_count += 1
                if haskey(parent_dict, child)
                    # ADD SPECIAL CASE HANDLING HERE
                else
                    parent_dict[child] = node_count
                end
                println("Node #$(length(new_nds)), (c#, a#, p#, n#, cn, an): $((child, active_node_num, parent_num, node_count, cn, added_node))")
            end
        end
        println("post queue: $(queue)")
    end
    return new_nds
end

function remove_subexpr_children!(expr::_NonlinearExprData)
    nd = expr.nd
    adj = adjmat(nd)
    children_arr = rowvals(adj)
    @inbounds first_node = nd[1]
    new_nds = NodeData[first_node]

    parent_dict = Dict{Int,Int}(1 => 1)
    queue = Tuple{Int,Int}[(1,-1)]
    node_count = 1

    while ~isempty(queue)
        (node_num, prior_prt) = popfirst!(queue)
        @inbounds active_node = nd[node_num]
        if (active_node.nodetype !== SUBEXPRESSION &&
            active_node.nodetype !== MOIVARIABLE &&
            active_node.nodetype !== VARIABLE &&
            active_node.nodetype !== VALUE)
            @inbounds children_idx = nzrange(adj, node_num)
            if (length(children_idx) > 0) # has any children
                for child in children_idx
                    @inbounds idx = children_arr[child]
                    @inbounds cn = nd[idx]
                    @inbounds parent_num = parent_dict[node_num]
                    push!(queue, (idx, node_num))
                    push!(new_nds, NodeData(cn.nodetype, cn.index, parent_num))
                    node_count += 1
                    parent_dict[idx] = node_count
                end
            end
        end
    end
    expr.nd = new_nds
end

function replace_subexpressions!(expr::_NonlinearExprData, mv_len::Int, n_expr0::Int)
    for i in 1:length(expr.nd)
        @inbounds nd = expr.nd[i]
        if nd.nodetype == CALL
            diff_indx = nd.index - USER_OPERATOR_ID_START
            if diff_indx >= 0
                expr.nd[i] = NodeData(SUBEXPRESSION, n_expr0 + diff_indx + 1, nd.parent)
            end
        elseif nd.nodetype == CALLUNIVAR
            diff_indx = nd.index - USER_UNIVAR_OPERATOR_ID_START
            if diff_indx >= 0
                expr.nd[i] = NodeData(SUBEXPRESSION, n_expr0 + mv_len + diff_indx + 1, nd.parent)
            end
        end
    end
end

function add_subexpr_from_tape!(tape::Tape, jnlp_data)
    nd = tape_to_list(tape)
    const_values = tape.const_values
    subexpr = _NonlinearExprData(nd, const_values)
    push!(jnlp_data.nlexpr, subexpr)
end

function udf_loader!(x::AbstractOptimizer)
    reform_flag = in(:presolve_flatten_flag, [fieldnames(typeof(x))...])
    if reform_flag
        reform_flag &= getfield(x, :presolve_flatten_flag)
    end

    evaluator = x._nlp_data.evaluator
    parameter_values = evaluator.parameter_values
    jnlp_data = evaluator.m.nlp_data
    user_registry = jnlp_data.user_operators
    user_reg_moe = user_registry.multivariate_operator_evaluator
    user_reg_uof = user_registry.univariate_operator_f
    n_expr0 = length(jnlp_data.nlexpr)

    # extracts tape from multivariate udf functions and creates subexpressions
    for ufe in user_reg_moe
        tape = trace_script(ufe.f, ufe.len)
        add_subexpr_from_tape!(tape, jnlp_data)
    end
    mv_len = length(user_reg_moe)
    for f in user_reg_uof
        tape = trace_script(f, 1)
        add_subexpr_from_tape!(tape, jnlp_data)
    end

    # replaces references in expressions to udfs with reference to subexpr
    # and remove any children of subexpressions (since subexpressions are terminal nodes)
    for i in 1:n_expr0
        @inbounds expr = jnlp_data.nlexpr[i]
        replace_subexpressions!(expr, mv_len, n_expr0)
        remove_subexpr_children!(expr)
        x.presolve_flatten_flag && flatten_expression!(expr, parameter_values)
    end
    if (jnlp_data.nlobj !== nothing)
        replace_subexpressions!(jnlp_data.nlobj, mv_len, n_expr0)
        remove_subexpr_children!(jnlp_data.nlobj)
        x.presolve_flatten_flag && flatten_expression!(jnlp_data.nlobj, parameter_values)
    end
    constr_len = length(jnlp_data.nlconstr)
    for i in 1:constr_len
        @inbounds constr = jnlp_data.nlconstr[i]
        replace_subexpressions!(constr.terms, mv_len, n_expr0)
        remove_subexpr_children!(constr.terms)
        x.presolve_flatten_flag && flatten_expression!(constr.terms, parameter_values)
    end

    # void previously defined udfs
    jnlp_data.user_operators = UserOperatorRegistry()
    jnlp_data.largest_user_input_dimension = 0
    evaluator.m.nlp_data = jnlp_data
    evaluator.eval_objective_timer = 0.0
    x._nlp_data = NLPBlockData(x._nlp_data.constraint_bounds, evaluator, x._nlp_data.has_objective)

    # reinitialize evaluator
    features = features_available(x._nlp_data.evaluator)
    init_feat = [:Grad, :Hess]
    num_nlp_constraints = length(x._nlp_data.constraint_bounds)
    num_nlp_constraints > 0 && push!(init_feat, :Jac)
    initialize(x._nlp_data.evaluator, init_feat)
    return
end
