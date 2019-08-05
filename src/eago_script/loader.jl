#=
scripts which take udfs, applies tracer, and converts them to JuMP nlexprssions
key components... transform forward DAG to reverse... then load model
takes... and reconstructs the nlp_data field appropriately
=#

# NOT A BREADFIRST SEARCH... EXHAUSTIVE SEARCH OF QUEUE MOSTLY DONE...
function tape_to_list(tape::Tape)

    len = length(tape.nd)
    last_node = tape.nd[len]
    new_nds = NodaData[NodeData(last_node.nodetype, last_node.index, -1)]

    queue = Tuple{Int,Int}[]
    push!(queue, (len, -1))
    while ~isempty(queue)
        (node_num, prior_prt) = popfirst!(queue)
        if (tape.nd[node_num].children[1] > 0)
            for child in tape.nd[node_num].children
                push!(queue, (child, node_num))
                cn = tape.nd[child]
                push!(new_nds, NodeData(cn.nodetype, cn.index, node_num))
            end
        end
    end
    return new_nds
end

function udf_loader!(m::MathOptInterface.AbstractOptimizer)

    nlp_block = m.nlp_data
    evaluator = nlp_block.evaluator
    user_registry = evaluator.user_operators
    jnlp_data = evaluator.m.nlp_data

    # extracts tape from multivariate udf functions and creates subexpressions
    for ufe in user_registry.multivariate_operator_evaluator
        tape = trace_script(ufe.f, ufe.len)
        add_subexpr_from_tape!(tape, jnlp_data)
    end
    mv_len = length(user_registry.multivariate_operator_evaluator)

    # extracts tape from univariate udf functions and creates subexpressions
    for f in user_registry.univariate_operator_f
        tape = trace_script(f, 1)
        add_subexpr_from_tape!(tape, jnlp_data)
    end

    # replaces references in expressions to udfs with reference to subexpr
    # (assumes no nested udfs)
    for expr in jnlp_data.nlexpr[1:n_expr0]
        for i in 1:length(expr.nd)
            if expr.nd[i].nodetype == CALL
                diff_indx = nd[i].index - USER_OPERATOR_ID_START
                if diff_indx >= 0
                    expr.nd[i] = NodeData(SUBEXPRESSION, n_expr0 + diff_indx + 1, expr.nd[i].parent)
                end
            elseif nd.nodetype == CALLUNIVAR
                diff_indx = nd[i].index - USER_UNIVAR_OPERATOR_ID_START
                if diff_indx >= 0
                    expr.nd[i] = NodeData(SUBEXPRESSION, n_expr0 + mv_len + diff_indx + 1, expr.nd[i].parent)
                end
            end
        end
    end

    # replaces references in objective to udfs with reference to subexpr
    for i in 1:length(jnlp_data.nlobj.nd)
        if jnlp_data.nlobj.nd[i].nodetype == CALL
            diff_indx = jnlp_data.nlobj.nd[i].index - USER_OPERATOR_ID_START
            if diff_indx >= 0
                jnlp_data.nlobj.nd[i] = NodeData(SUBEXPRESSION, n_expr0 + diff_indx + 1, jnlp_data.nlobj.nd[i].parent)
            end
        elseif nd.nodetype == CALLUNIVAR
            diff_indx = jnlp_data.nlobj.nd[i].index - USER_UNIVAR_OPERATOR_ID_START
            if diff_indx >= 0
                jnlp_data.nlobj.nd[i] = NodeData(SUBEXPRESSION, n_expr0 + mv_len + diff_indx + 1, jnlp_data.nlobj.nd[i].parent)
            end
        end
    end

    # replaces references in constraints to udfs with reference to subexpr
    for expr in jnlp_data.nlconstr
        for i in 1:length(expr.nd)
            if expr.nd[i].nodetype == CALL
                diff_indx = nd[i].index - USER_OPERATOR_ID_START
                if diff_indx >= 0
                    expr.nd[i] = NodeData(SUBEXPRESSION, n_expr0 + diff_indx + 1, expr.nd[i].parent)
                end
            elseif nd.nodetype == CALLUNIVAR
                diff_indx = nd[i].index - USER_UNIVAR_OPERATOR_ID_START
                if diff_indx >= 0
                    expr.nd[i] = NodeData(SUBEXPRESSION, n_expr0 + mv_len + diff_indx + 1, expr.nd[i].parent)
                end
            end
        end
    end

    # void previously defined udfs
    m.nlp_data.evaluator.user_operators = UserOperatorRegistry()

    # reinitialize evaluator
    features = MOI.features_available(evaluator)
    has_hessian = (:Hess in features)
    init_feat = [:Grad, :ExprGraph]
    num_nlp_constraints > 0 && push!(init_feat, :Jac)
    MOI.initialize(evaluator, init_feat)
end

#=
function tape_to_lists(tape::Tape; new_exprs = _NonlinearExprData[])
    const_values = tape.const_values
    len = length(nd)
    nd = Vector{NodeData}(undef, len)
    parent_dict = Dict{Int,Int}()
    new_tape_dict = Dict{Int,Tape}()
    shared_dict = Dict{Int,Int}() # original child to new expression
    shared_child_count = 1
    for i in len:-1:1
        # assign terminal value parent == -1
        if i == len
           nd[i] = NodeData(nd.nodetype, nd.index, -1)
        end
        for child in nd.nodetype.children
            if haskey(parent_dict, child)
                if haskey(new_tape_dict, child)
                error("node has two parents... ")
                # TODO: child has two parents which JuMP does't allow ->
                # make child node a subexpression and have both parents reference it
                # nd, const_values, new_exprs
            else
                parent_dict[child] = i
            end
        end
        parent_loc = parent_dict[i]
        nd[i] = NodeData(nd.nodetype, nd.index, parent_loc)
    end

    # reverses node array adjusting parent nodes appropriately
    nd_out = Vector{NodeData}(undef, len)
    count = 1
    for i=length(nd):-1:1
        if (parent_loc[i] > 0)
            nd_out[count] = NodeData(nd.nodetype[i], nd.index[i], length(nd) + i - parent_loc[i])
        else
            nd_out[count] = NodeData(nd.nodetype[i], nd.index[i], parent_loc[i])
        end
        count += 1
    end

    return nd, const_values, new_exprs
end

function add_subexpr_from_tape!(tape::Tape, jnlp_data)
    nd, const_values, new_exprs = tape_to_lists(tape)
    subexpr = JuMP._NonlinearExprData(nd, const_values)
    append!(jnlp_data, new_exprs)
    push!(jnlp_data, subexpr)
end

function udf_loader!(m::MathOptInterface.AbstractOptimizer)

    nlp_block = m.nlp_data   # TODO: search optimizer for nlp_data field type... pull data from that if it exists
    evaluator = nlp_block.evaluator
    user_registry = evaluator.user_operators
    jnlp_data = evaluator.m.nlp_data

    n_expr0 = length(jnlp_data.nlexpr)

    # extracts tape from multivariate udf functions and creates subexpressions
    for ufe in user_registry.multivariate_operator_evaluator
        tape = trace_script(ufe.f, ufe.len)
        add_subexpr_from_tape!(tape, jnlp_data)
    end
    mv_len = length(user_registry.multivariate_operator_evaluator)

    # extracts tape from univariate udf functions and creates subexpressions
    for f in user_registry.univariate_operator_f
        tape = trace_script(f, 1)
        add_subexpr_from_tape!(tape, jnlp_data)
    end

    # replaces references in expressions to udfs with reference to subexpr
    # (assumes no nested udfs)
    for expr in jnlp_data.nlexpr[1:n_expr0]
        for i in 1:length(expr.nd)
            if expr.nd[i].nodetype == CALL
                diff_indx = nd[i].index - USER_OPERATOR_ID_START
                if diff_indx >= 0
                    expr.nd[i] = NodeData(SUBEXPRESSION, n_expr0 + diff_indx + 1, expr.nd[i].parent)
                end
            elseif nd.nodetype == CALLUNIVAR
                diff_indx = nd[i].index - USER_UNIVAR_OPERATOR_ID_START
                if diff_indx >= 0
                    expr.nd[i] = NodeData(SUBEXPRESSION, n_expr0 + mv_len + diff_indx + 1, expr.nd[i].parent)
                end
            end
        end
    end

    # replaces references in objective to udfs with reference to subexpr
    for i in 1:length(jnlp_data.nlobj.nd)
        if jnlp_data.nlobj.nd[i].nodetype == CALL
            diff_indx = jnlp_data.nlobj.nd[i].index - USER_OPERATOR_ID_START
            if diff_indx >= 0
                jnlp_data.nlobj.nd[i] = NodeData(SUBEXPRESSION, n_expr0 + diff_indx + 1, jnlp_data.nlobj.nd[i].parent)
            end
        elseif nd.nodetype == CALLUNIVAR
            diff_indx = jnlp_data.nlobj.nd[i].index - USER_UNIVAR_OPERATOR_ID_START
            if diff_indx >= 0
                jnlp_data.nlobj.nd[i] = NodeData(SUBEXPRESSION, n_expr0 + mv_len + diff_indx + 1, jnlp_data.nlobj.nd[i].parent)
            end
        end
    end

    # replaces references in constraints to udfs with reference to subexpr
    for expr in jnlp_data.nlconstr
        for i in 1:length(expr.nd)
            if expr.nd[i].nodetype == CALL
                diff_indx = nd[i].index - USER_OPERATOR_ID_START
                if diff_indx >= 0
                    expr.nd[i] = NodeData(SUBEXPRESSION, n_expr0 + diff_indx + 1, expr.nd[i].parent)
                end
            elseif nd.nodetype == CALLUNIVAR
                diff_indx = nd[i].index - USER_UNIVAR_OPERATOR_ID_START
                if diff_indx >= 0
                    expr.nd[i] = NodeData(SUBEXPRESSION, n_expr0 + mv_len + diff_indx + 1, expr.nd[i].parent)
                end
            end
        end
    end

    # void previously defined udfs
    m.nlp_data.evaluator.user_operators = UserOperatorRegistry()

    # reinitialize evaluator
    features = MOI.features_available(evaluator)
    has_hessian = (:Hess in features)
    init_feat = [:Grad, :ExprGraph]
    num_nlp_constraints > 0 && push!(init_feat, :Jac)
    MOI.initialize(evaluator, init_feat)
end


# NodeData(SUBEXPRESSION, loc, parent)
# DFS for colored nodes (DONE)
# Add all nodes at each colored node to working tape
# delete from prior tape if nodes has a parents outside new tape and add to new tape)
#

function adjmat(x::Tape)
    len = length(x.nd)
    adj = spzeros(Int, len.nd, len.nd)
    for i in len:-1:1
        for child in x.nd[i].children
            adj[child,i] = 1
        end
    end
    return adj
end

function tape_to_JuMP(x::Tape)
    # colors nodes by number of parents
    adj = adjmat(x)
    parent_number_dict = Dict{Int,Int}([i => sum(a[i,:]) for i in 1:length(x.nd)])

    # create storage between
    num_new_exprs = sum(x -> x[2] > 1, parent_number_dict)
    new_exprs = _NonlinearExprData[_NonlinearExprData(NodeData[], Float64[]) for i in 1:num_new_exprs]
    new_exprs_location = keys(filter(x -> x[2] > 1, parent_number_dict))

    labelled_nodes alue = nd::Vector{NodeInfo}

    # sets up storage for dfs
    dfs_core()
    # performs dfs
end

function dfs_strategy(x::Tape, 1)
    tape_dict = Dict{Int,NodeInfo}([i => x.nd[i] for i in 1:length(x.nd)])
    visited_dict = Dict{Int,Bool}([false for i in 1:length(x.nd)])

    # locate multiparents via depth first search
    for i

end
=#
