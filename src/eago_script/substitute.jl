"""
    Template_Node

A structure which holds a symbol indicating whether the node is an operator,
a number, or an expression `type`, a value which identifies the function or
symbol `value`, potentially a numeric value `num_value`, and a check that can
be run to verify the node is correct `check`.
"""
struct Template_Node <: Any
    type::Symbol # op, num, expr
    value::Symbol
    num_value::Float64
    check::Function
end

"""
    Template_Graph

Holds a list of Template_Nodes, set of directed edges, lengths, an adjacency
matrix and the number of children.
"""
struct Template_Graph <: Any
    nd::Vector{Template_Node}
    dag::Vector{Pair{Int,Int}}
    ndlen::Int
    daglen::Int
    adj::SparseMatrixCSC{Bool,Int}
    num_children::Vector{Int}
end
function Template_Graph(nd::Dict{Int,Template_Node}, dag::Vector{Pair{Int,Int}})
    ndlen = length(nd)
    daglen = length(dag)
    ndlist = [nd[i] for i in 1:ndlen]
    adj = spzeros(Bool, ndlen, ndlen)
    for i in 1:daglen
        @inbounds p = dag[i]
        @inbounds x = p[1]
        @inbounds y = p[2]
        @inbounds adj[x, y] = true
    end
    num_children = Int[length(nzrange(adj, i)) for i in 1:ndlen]
    Template_Graph(ndlist, dag, ndlen, daglen, adj, num_children)
end

always_true(x) = true
function Template_Node(type::Symbol, value::Symbol; check::Function = always_true)
    @assert (type == :op || type == :expr || type == :num)
    Template_Node(type, value, NaN, check)
end
function Template_Node(type::Symbol, num_value::Float64; check::Function = always_true)
    @assert (type == :num)
    Template_Node(type, :null, num_value, check)
end

const DAG_PATTERNS = Template_Graph[]
const DAG_SUBSTITUTIONS = Template_Graph[]
const DAG_SPDICT = Dict{Int,Int}[]
const DAG_LENGTHS = Int[0,0]

#=
conventions for substition, the expression to be checked always appears at key 1
in the Template_Graph and operations are ordered from low value to high value left to right
so if 1 is a -, and 4 => 1, 3 => 1 then the expression is 4 - 3
=#
"""
    register_substitution!

Specifies that the `src::Template_Graph` should be subsituted out for the
`trg::Template_Graph`.

Conventions for substition, the expression to be checked always appears at key 1
in the Template_Graph and operations are ordered from low value to high value left to right
so if 1 is a -, and 4 => 1, 3 => 1 then the expression is 4 - 3
"""
function register_substitution!(src::Template_Graph, trg::Template_Graph)
    push!(DAG_PATTERNS, src)
    push!(DAG_SUBSTITUTIONS, trg)
    DAG_LENGTHS[1] += 1
    DAG_LENGTHS[2] += 1
    d = Dict{Int,Int}()
    for (i, src_nd) in enumerate(src.nd)
        src_nd_type = src_nd.type
        if (src_nd_type === :expr || src_nd_type === :num)
            for (j, trg_nd) in enumerate(trg.nd)
                if trg_nd.type === src_nd_type
                    if src_nd.value === trg_nd.value
                        @inbounds d[j] = i
                        break
                    end
                end
            end
        end
    end
    push!(DAG_SPDICT, d)
end

function matching_info(x::Template_Node, y::NodeData,
                       const_values::Vector{Float64}, parameter_values::Vector{Float64})
    match_flag = false
    if (x.type == :op)
        if (y.nodetype == CALL)
            if haskey(operator_to_id, x.value)
                @inbounds indx_id = operator_to_id[x.value]
                if indx_id == y.index
                    match_flag = true
                end
            else
                match_flag = false
            end
        elseif (y.nodetype == CALLUNIVAR)
            if haskey(univariate_operator_to_id, x.value)
                @inbounds indx_id = univariate_operator_to_id[x.value]
                if indx_id == y.index
                    match_flag = true
                end
            else
                match_flag = false
            end
        end
    elseif (x.type == :num)
        if (y.nodetype == VALUE)
            @inbounds num_value = const_values[y.index]
        elseif (y.nodetype == PARAMETER)
            @inbounds num_value = parameter_values[y.index]
        else
            match_flag = false
        end
        if (y.nodetype == VALUE) || (y.nodetype == PARAMETER)
            if x.check(num_value)
                match_flag = true
            else
                match_flag = false
            end
        end
    elseif (x.type == :expr)
        match_flag = true
    end
    return match_flag
end

function is_match(pattern::Template_Graph, indx::Int, nd::Vector{NodeData}, dag_adj::SparseMatrixCSC{Bool,Int},
                  const_values::Vector{Float64}, parameter_values::Vector{Float64})
    match_flag = true
    match_dict = Dict{Int,Int}()
    pattern_length = pattern.ndlen
    dag_length = pattern.daglen
    pattern_adj = pattern.adj
    pat_children_arr = rowvals(pattern_adj)
    dag_children_arr = rowvals(dag_adj)

    # do a breadth first search with paired template, nd data,
    # if any pair of children fail then
    pindx_initial = 1
    queue = Tuple{Int,Int}[(pindx_initial, indx)]
    while (~isempty(queue) && (match_flag == true))
        (num_pat, num_dag) = popfirst!(queue)
        @inbounds patt_nd = pattern.nd[num_pat]
        @inbounds dag_nd = nd[num_dag]
        match_info_flag = matching_info(patt_nd, dag_nd, const_values, parameter_values)
        if match_info_flag
            if patt_nd.type === :expr
                @inbounds match_dict[num_pat] = num_dag
            elseif (patt_nd.type === :num) && (patt_nd.num_value === NaN)
                @inbounds match_dict[num_pat] = num_dag
            end
            @inbounds pat_children_idx = nzrange(pattern_adj, num_pat)
            pat_length = length(pat_children_idx)
            if pat_length > 0
                @inbounds dag_children_idx = nzrange(dag_adj, num_dag)
                if pat_length == length(dag_children_idx)
                    for i in 1:pat_length
                        @inbounds dchild = dag_children_idx[i]
                        @inbounds pchild = pat_children_idx[i]
                        @inbounds didx = dag_children_arr[dchild]
                        @inbounds pidx = pat_children_arr[pchild]
                        push!(queue, (pidx, didx))
                    end
                else
                    match_flag = false
                    break
                end
            end
        else
            match_flag = false
            break
        end
    end
    return match_flag, match_dict
end

function find_match(indx::Int, nd::Vector{NodeData}, adj::SparseMatrixCSC{Bool,Int},
                    const_values::Vector{Float64}, parameter_values::Vector{Float64})

    flag = false
    pattern_number = -1
    match_dict = Dict{Int,Int}()
    @inbounds sub_len = DAG_LENGTHS[1]
    for i in 1:sub_len
        @inbounds pattern = DAG_PATTERNS[i]
        inner_flag, match_dict = is_match(pattern, indx, nd, adj,
                                          const_values, parameter_values)
        if inner_flag
            flag = true
            pattern_number = i
            break
        end
    end
    return flag, pattern_number, match_dict
end

#=
Takes a template node and makes the appropriate JuMP node, takes the parent index,
number of child for a pattern element, constant storage vector and it's length
=#
function op_node_to_dag!(x::Template_Node, parent::Int, child_len::Int)
    if child_len > 1
        @inbounds op = operator_to_id[x.value]
        node = NodeData(CALL, op, parent)
    else
        @inbounds op = univariate_operator_to_id[x.value]
        node = NodeData(CALLUNIVAR, op, parent)
    end
    return node
end

function bfs_expr_add!(new_nds::Vector{NodeData}, node_count::Int, num_prt::Int,
                       parent_dict::Dict{Int,Int}, match_dict::Dict{Int,Int},
                       expr_loc::Int, nd::Vector{NodeData}, adj::SparseMatrixCSC{Bool,Int},
                       children_arr::Vector{Int})
    queue = Tuple{Int,Int}[(expr_loc, num_prt)]
    inner_node_count = node_count
    while ~isempty(queue)
        (node_num, prior_prt) = popfirst!(queue) # pop node
        @inbounds active_node = nd[node_num] # store node
        new_node = NodeData(active_node.nodetype, active_node.index, prior_prt)
        inner_node_count += 1 # update node count
        @inbounds parent_dict[num_prt] = inner_node_count
        push!(new_nds, new_node)
        if (active_node.nodetype !== SUBEXPRESSION &&
            active_node.nodetype !== MOIVARIABLE &&
            active_node.nodetype !== VARIABLE &&
            active_node.nodetype !== VALUE)
            @inbounds children_idx = nzrange(adj, node_num)
            if (length(children_idx) > 0)
                for child in children_idx
                    @inbounds idx = children_arr[child]
                    push!(queue, (idx, inner_node_count))
                end
            end
        end
    end
    inner_node_count
end

# we assume a tree structure, so if we don't load child nodes,
# then they are effectively deleted
function substitute!(match_num::Int, node_num::Int, prior_prt::Int, nd::Vector{NodeData},
                     const_list::Vector{Float64}, const_len::Int, node_count::Int,
                     parent_dict::Dict{Int,Int}, match_dict::Dict{Int,Int},
                     queue::Vector{Tuple{Int,Int}}, new_nds::Vector{NodeData},
                     adj::SparseMatrixCSC{Bool,Int}, children_arr::Vector{Int})
    @inbounds subs_template = DAG_SUBSTITUTIONS[match_num]
    @inbounds subs_patt_dict = DAG_SPDICT[match_num]
    subs_adj = subs_template.adj
    subs_children_arr = rowvals(subs_adj)
    queue = Tuple{Int,Int,Int}[(prior_prt, 1, -1)] # node_num, prior_prt, dag_num, dag_prt
    inner_node_count = node_count
    while ~isempty(queue)
        (num_prt, num_sub, num_sub_prt) = popfirst!(queue)
        @inbounds active_node = subs_template.nd[num_sub]
        active_type = active_node.type
        if active_type === :op
            @inbounds active_node_children = subs_template.num_children[num_sub]
            node = op_node_to_dag!(active_node, num_prt, active_node_children)
            push!(new_nds, node)
            inner_node_count += 1
            @inbounds parent_dict[num_prt] = inner_node_count
            @inbounds subs_children_idx = nzrange(subs_adj, num_sub)
            if ~isempty(subs_children_idx)
                for child in subs_children_idx
                    @inbounds sub_child = subs_children_arr[child]
                    push!(queue, (inner_node_count, sub_child, num_sub))
                end
            end
        elseif active_type === :num
            if (active_node.num_value !== NaN)
                push!(const_list, active_node.num_value)
                const_len += 1
                inner_node_count += 1
                @inbounds parent_dict[num_prt] = inner_node_count
                node = NodeData(VALUE, const_len, num_prt)
                push!(new_nds, node)
            else
                @inbounds pindx = subs_patt_dict[num_sub]
                @inbounds expr_loc = match_dict[pindx]
                inner_node_count += 1
                @inbounds parent_dict[num_prt] = inner_node_count
                @inbounds node = NodeData(VALUE, nd[expr_loc].index, num_prt)
                push!(new_nds, node)
            end
        elseif active_type === :expr # Need to only use
            @inbounds pindx = subs_patt_dict[num_sub]
            @inbounds expr_loc = match_dict[pindx]
            temp_node_count = bfs_expr_add!(new_nds, inner_node_count, num_prt,
                                            parent_dict, match_dict, expr_loc,
                                            nd, adj, children_arr)
            inner_node_count = temp_node_count
        end
    end
    return inner_node_count, const_len
end

# searchs through expression breadth first search that short-cirucits
"""
    flatten_expression!

Flattens (usually) the dag by making all registered substitutions for the
expression `expr::_NonlinearExprData`.
"""
function flatten_expression!(expr::_NonlinearExprData, parameter_values::Vector{Float64})
    nd = expr.nd
    adj = adjmat(nd)
    children_arr = rowvals(adj)
    node_count = 0
    const_list = expr.const_values
    const_len = length(const_list)
    parent_dict = Dict{Int,Int}(-1 => -1)
    queue = Tuple{Int,Int}[(1,-1)]
    new_nds = NodeData[]
    while ~isempty(queue)
        (node_num, prior_prt) = popfirst!(queue)
        @inbounds active_node = nd[node_num]
        if (active_node.nodetype !== SUBEXPRESSION &&
            active_node.nodetype !== MOIVARIABLE &&
            active_node.nodetype !== VARIABLE &&
            active_node.nodetype !== VALUE)
            is_match, match_num, match_dict = find_match(node_num, nd, adj, const_list, parameter_values)
            if ~is_match
                @inbounds parent_num = parent_dict[prior_prt]
                push!(new_nds, NodeData(active_node.nodetype, active_node.index, parent_num))
                node_count += 1
                @inbounds parent_dict[node_num] = node_count
                @inbounds children_idx = nzrange(adj, node_num) # ADD CHILDREN
                if (length(children_idx) > 0)
                    for child in children_idx
                        @inbounds idx = children_arr[child]
                        push!(queue, (idx, node_num))
                    end
                end
             else
                node_count, const_len = substitute!(match_num, node_num, prior_prt,
                                                    nd, const_list, const_len,
                                                    node_count, parent_dict, match_dict,
                                                    queue, new_nds, adj, children_arr)
            end
        else
            push!(new_nds, NodeData(active_node.nodetype, active_node.index, prior_prt))
            node_count += 1
            @inbounds parent_dict[node_num] = node_count
        end
    end
    expr.nd = new_nds
end

"""
    dag_flattening!

Flattens (usually) the dag by making all registered substitutions for every
nonlinear term in the Optimizer.
"""
function dag_flattening!(x::T) where T <: AbstractOptimizer
    nlp_data = x._nlp_data.evaluator.m.nlp_data
    params = nlp_data.nlparamvalues
    if ~isnothing(nlp_data.nlobj)
        flatten_expression!(nlp_data.nlobj, params)
    end
    for i in 1:length(nlp_data.nlconstr)
        flatten_expression!(nlp_data.nlconstr[i].terms, params)
    end
    for i in 1:length(nlp_data.nlexpr)
        flatten_expression!(nlp_data.nlexpr[i], params)
    end
    return
end
