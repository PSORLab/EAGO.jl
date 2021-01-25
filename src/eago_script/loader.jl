# Copyright (c) 2018: Matthew Wilhelm & Matthew Stuber.
# This code is licensed under MIT license (see LICENSE.md for full details)
#############################################################################
# EAGO
# A development environment for robust and global optimization
# See https://github.com/PSORLab/EAGO.jl
#############################################################################
# src/eago_script/loader.jl
# Utilities used to load user-defined functions into JuMP tapes.
#############################################################################

"""
"""
function tape_to_list(tape::Tape)
    len = length(tape.nd)
    @inbounds last_node = tape.nd[len]
    new_nds = NodeData[NodeData(last_node.nodetype, last_node.index, -1)]

    queue = Tuple{Int,Int}[(len, -1)]
    parent_dict = Dict{Int,Int}(len => 1) # starting node is 1
    node_count = 1

    while !isempty(queue)
        (active_node_num, prior_prt) = popfirst!(queue)
        @inbounds active_node = tape.nd[active_node_num]
        @inbounds active_node_child1 = active_node.children[1]

        if active_node_child1 > 0
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
            end
        end
    end

    return new_nds
end

"""
"""
function remove_subexpr_children!(expr::_NonlinearExprData)
    nd = expr.nd
    adj = adjmat(nd)
    children_arr = rowvals(adj)
    @inbounds first_node = nd[1]
    new_nds = NodeData[first_node]

    parent_dict = Dict{Int,Int}(1 => 1)
    queue = Tuple{Int,Int}[(1,-1)]
    node_count = 1

    while !isempty(queue)
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

    return nothing
end

"""
"""
function replace_subexpressions!(expr::_NonlinearExprData, mv_len::Int, nlexpr_count::Int)
    for i = 1:length(expr.nd)
        nd = @inbounds expr.nd[i]
        if nd.nodetype === CALL
            diff_indx = nd.index - USER_OPERATOR_ID_START
            if diff_indx >= 0
                expr.nd[i] = NodeData(SUBEXPRESSION, nlexpr_count + diff_indx + 1, nd.parent)
            end
        elseif nd.nodetype === CALLUNIVAR
            diff_indx = nd.index - USER_UNIVAR_OPERATOR_ID_START
            if diff_indx >= 0
                expr.nd[i] = NodeData(SUBEXPRESSION, nlexpr_count + mv_len + diff_indx + 1, nd.parent)
            end
        end
    end

    return nothing
end

"""
"""
function add_subexpr_from_tape!(tape::Tape, jnlp_data)
    nd = tape_to_list(tape)
    const_values = tape.const_values
    subexpr = _NonlinearExprData(nd, const_values)
    push!(jnlp_data.nlexpr, subexpr)

    return nothing
end

"""
"""
function udf_loader!(x::AbstractOptimizer)
    reform_flag = in(:presolve_flatten_flag, [fieldnames(typeof(x))...])
    if reform_flag
        reform_flag &= getfield(x, :presolve_flatten_flag)
    end

    evaluator = x._working_problem._nlp_data.evaluator
    parameter_values = evaluator.parameter_values
    nlp_data = evaluator.m.nlp_data
    user_registry = nlp_data.user_operators

    # extracts tape from multivariate udf functions and creates subexpressions
    multi_op_eval = user_registry.multivariate_operator_evaluator
    for mul_eval in multi_op_eval
        tape = trace_script(mul_eval.f, mul_eval.len)
        add_subexpr_from_tape!(tape, nlp_data)
    end
    mv_len = length(multi_op_eval)

    uni_op_f = user_registry.univariate_operator_f
    for f in uni_op_f
        tape = trace_script(f, 1)
        add_subexpr_from_tape!(tape, nlp_data)
    end

    # replaces references in expressions to udfs with reference to subexpr
    # and remove any children of subexpressions (since subexpressions are terminal nodes)
    nlexpr = nlp_data.nlexpr
    nlexpr_count = length(nlexpr)
    for i = 1:nlexpr_count
        expr = @inbounds nlp_data.nlexpr[i]
        replace_subexpressions!(expr, mv_len, nlexpr_count)
        remove_subexpr_children!(expr)
        x.presolve_flatten_flag && flatten_expression!(expr, parameter_values)
    end

    if nlp_data.nlobj !== nothing
        replace_subexpressions!(nlp_data.nlobj, mv_len, nlexpr_count)
        remove_subexpr_children!(nlp_data.nlobj)
        x.presolve_flatten_flag && flatten_expression!(nlp_data.nlobj, parameter_values)
    end

    nlconstr = nlp_data.nlconstr
    constr_len = length(nlconstr)
    for i = 1:constr_len
        constr = @inbounds nlconstr[i]
        replace_subexpressions!(constr.terms, mv_len, nlexpr_count)
        remove_subexpr_children!(constr.terms)
        x.presolve_flatten_flag && flatten_expression!(constr.terms, parameter_values)
    end

    # void previously defined udfs
    nlp_data.user_operators = UserOperatorRegistry()
    nlp_data.largest_user_input_dimension = 0
    evaluator.m.nlp_data = nlp_data
    evaluator.eval_objective_timer = 0.0
    x._nlp_data = NLPBlockData(x._nlp_data.constraint_bounds, evaluator, x._nlp_data.has_objective)

    # reinitialize evaluator
    features = features_available(x._nlp_data.evaluator)
    init_feat = Symbol[:Grad, :Hess]
    num_nlp_constraints = length(x._nlp_data.constraint_bounds)
    num_nlp_constraints > 0 && push!(init_feat, :Jac)
    MOI.initialize(x._nlp_data.evaluator, init_feat)

    return nothing
end
