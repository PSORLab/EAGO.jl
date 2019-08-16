struct Template_Node <: Any
    type::Symbol # op, num, expr
    value::Symbol
    num_value::Float64
    check::Function
end
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
        adj[x, y] = true
    end
    num_children = Int[length(nzrange(adj, i)) for i in 1:ndlen]
    # compute number of children per node
    #println("ndlist: $ndlist")
#    println("dag: $dag")
    #println("ndlen: $ndlen")
    #println("daglen: $daglen")
    #println("adj: $adj")
    #println("num_children: $num_children")
    Template_Graph(ndlist, dag, ndlen, daglen, adj, num_children)
end

always_true(x) = true
function Template_Node(type::Symbol, value::Symbol; check::Function = always_true)
    @assert (type == :op || type == :expr || type == :num)
    Template_Node(type, value, 0.0, check)
end
function Template_Node(type::Symbol, num_value::Float64; check::Function = always_true)
    @assert (type == :num)
    Template_Node(type, :null, num_value, check)
end

const DAG_PATTERNS = Template_Graph[]
const DAG_SUBSTITUTIONS = Template_Graph[]
const DAG_LENGTHS = Int[0,0]

#=
conventions for substition, the expression to be checked always appears at key 1
in the Template_Graph and operations are ordered from low value to high value left to right
so if 1 is a -, and 4 => 1, 3 => 1 then the expression is 4 - 3
=#
function register_substitution!(src::Template_Graph, trg::Template_Graph)
    push!(DAG_PATTERNS, src)
    push!(DAG_SUBSTITUTIONS, trg)
    DAG_LENGTHS[1] += 1
    DAG_LENGTHS[2] += 1
end

# SHOULD BE DONE
function matching_info(x::Template_Node, y::NodeData,
                       const_values::Vector{Float64}, parameter_values::Vector{Float64})
    match_flag = false
    if (x.type == :op)
        if (y.nodetype == CALL)
            if haskey(operator_to_id, x.value)
                if (operator_to_id[x.value] == y.index)
                    match_flag = true
                end
            else
                match_flag = false
            end
        elseif (y.nodetype == CALLUNIVAR)
            if haskey(univariate_operator_to_id, x.value)
                if univariate_operator_to_id[x.value] == y.index
                    match_flag = true
                end
            else
                match_flag = false
            end
        #elseif y.nodetype == SUBEXPRESSION TODO: Add subexpression handling latter
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
    else
        error("invalid type symbol")
    end
    return match_flag
end

# SHOULD BE DONE with the exception of matches that have no shared numbers
# and/or expression symbols (DEFINITELY AN ISSUE RIGHT NOW)
function is_match(pattern::Template_Graph, indx::Int, nd::Vector{NodeData}, dag_adj::SparseMatrixCSC{Bool,Int},
                  const_values::Vector{Float64}, parameter_values::Vector{Float64})
    #println("*** IS MATCH START ***")
    match_flag = true
    match_dict = Dict{Int,Int}()
    # make pattern adjacency matrix
    pattern_length = pattern.ndlen
    dag_length = pattern.daglen
    pattern_adj = pattern.adj
    pat_children_arr = rowvals(pattern_adj)
    dag_children_arr = rowvals(dag_adj)

    # do a breadth first search with paired template, nd data,
    # if any pair of children fail then
    pindx_initial = 1
    queue = Tuple{Int,Int}[(pindx_initial, indx)]
    depth = 0
    depth_max = 30
    while (~isempty(queue) && (match_flag == true) && (depth < depth_max))
        depth += 1
    #    println("depth: $depth")
    #    println("pre-pop queue: $queue")
        (num_pat, num_dag) = popfirst!(queue)
    #    println("num_pat: $num_pat")
    #    println("num_dag: $num_dag")
        @inbounds patt_nd = pattern.nd[num_pat]
        @inbounds dag_nd = nd[num_dag]
        #println("patt_nd: $patt_nd")
        #println("dag_nd: $dag_nd")
        match_info_flag = matching_info(patt_nd, dag_nd, const_values, parameter_values)
    #    println("match_info_flag: $match_info_flag")
        if match_info_flag
            #println("patt_nd.type: $(patt_nd.type)")
            if patt_nd.type == :expr
                match_dict[num_pat] = num_dag
            end
            @inbounds pat_children_idx = nzrange(pattern_adj, num_pat)
            pat_length = length(pat_children_idx)
            #println("pat_children_idx: $pat_children_idx")
            #println("pat_length: $pat_length")
            if pat_length > 0
                @inbounds dag_children_idx = nzrange(dag_adj, num_dag)
                #println("dag_children_idx: $dag_children_idx")
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
    #println("is match return: ")
    #println("match_flag: $match_flag")
    #println("match_dict: $match_dict")
    #println("*** IS MATCH END ***")
    return match_flag, match_dict
end

function find_match(indx::Int, nd::Vector{NodeData}, adj::SparseMatrixCSC{Bool,Int},
                    const_values::Vector{Float64}, parameter_values::Vector{Float64})

    #println("---FIND MATCH RETURN (START)---")
    flag = false
    pattern_number = -1
    match_dict = Dict{Int,Int}()
    @inbounds sub_len = DAG_LENGTHS[1]
    for i in 1:sub_len
        @inbounds pattern = DAG_PATTERNS[i]
        inner_flag, match_dict = is_match(pattern, indx, nd, adj, const_values, parameter_values)
        #println("inner_flag: $inner_flag")
        #println("match_dict: $match_dict")
        if inner_flag
            flag = true
            pattern_number = i
            break
        end
    end
    #println("flag $flag")
    #println("pattern_number $pattern_number")
    #println("match_dict $match_dict")
    #println("---FIND MATCH RETURN (END)---")
    return flag, pattern_number, match_dict
end

#=
Takes a template node and makes the appropriate JuMP node, takes the parent index,
number of child for a pattern element, constant storage vector and it's length
=#
function op_node_to_dag!(x::Template_Node, parent::Int, child_len::Int)
    if child_len > 1
        op = operator_to_id[x.value]
        node = NodeData(CALL, op, parent)
    else
        op = univariate_operator_to_id[x.value]
        node = NodeData(CALLUNIVAR, op, parent)
    end
    return node
end


# we assume a tree structure, so if we don't load child nodes,
# then they are effectively deleted
function substitute!(match_num::Int, node_num::Int, prior_prt::Int, nd::Vector{NodeData},
                     const_list::Vector{Float64}, const_len::Int, node_count::Int,
                     parent_dict::Dict{Int,Int}, match_dict::Dict{Int,Int}, queue::Vector{Tuple{Int,Int}},
                     new_nds::Vector{NodeData})
    subs_template = DAG_SUBSTITUTIONS[match_num]
    subs_adj = subs_template.adj
    subs_children_arr = rowvals(subs_adj)
    queue = Tuple{Int,Int,Int}[(prior_prt, 1, -1)] # node_num, prior_prt, dag_num, dag_prt
    depth = 0
    depth_max = 30
    println("begin substitute")
    println("match dict: $(match_dict)")
    println("starting new_nds: $new_nds")
    while ~isempty(queue) && (depth < depth_max)
        depth += 1
        println("substitute iteration #: $depth, length(queue): $(length(queue))")
        println("current new_nds: $new_nds")
        (num_prt, num_sub, num_sub_prt) = popfirst!(queue)
        active_node = subs_template.nd[num_sub]
        active_type = active_node.type
        println("active_type: $(active_type)")
        if active_type === :op
            active_node_children = subs_template.num_children[num_sub]
            println("active node childre: $(active_node_children)")
            node = op_node_to_dag!(active_node, num_prt, active_node_children)
            push!(new_nds, node)
            node_count += 1
            parent_dict[num_prt] = node_count
            @inbounds subs_children_idx = nzrange(subs_adj, num_sub)
            println("subs_children_idx: $(subs_children_idx)")
            if ~isempty(subs_children_idx)
                for child in subs_children_idx
                    sub_child = subs_children_arr[child]
                    println("sub_child = $sub_child")
                    push!(queue, (node_count, sub_child, num_sub)) # THE ERROR IS HERE....
                end
            end
        elseif active_type === :num
            push!(const_list, active_node.num_value)
            const_len += 1
            node_count += 1
            parent_dict[num_prt] = node_count
            node = NodeData(VALUE, const_len, num_prt)
            push!(new_nds, node)
        elseif active_type === :expr # Need to only use
            println("num_sub: $num_sub")
            expr_loc = match_dict[num_sub]
            println("expr_loc: $expr_loc")
            node_prior = nd[expr_loc]
            println("node_prior: $node_prior")
            node = NodeData(node_prior.nodetype, node_prior.index, num_prt)
            node_count += 1
            parent_dict[num_prt] = node_count
            push!(new_nds, node)
        end
    end
    println("end substitute")
    return node_count, const_len
end

# searchs through expression breadth first search that short-cirucits
# TODO: NEED TO CHECK, SHOULD BE DONE (how to shift )... no need to shift nodes
# since parents are reference from piror list
function flatten_expression!(expr::_NonlinearExprData, parameter_values::Vector{Float64})
    nd = expr.nd
    adj = adjmat(nd)
    children_arr = rowvals(adj)
    node_count = 1
    const_list = expr.const_values
    const_len = length(const_list)
    parent_dict = Dict{Int,Int}()
    queue = Tuple{Int,Int}[(1,-1)]
    new_nds = NodeData[]

    println("(flatten iteration) end depth = 0 check")
    depth = 0
    depth_max = 30
    while ~isempty(queue) && (depth < depth_max)
        depth += 1
        println("(flatten iteration) start depth = $depth check")
        (node_num, prior_prt) = popfirst!(queue)
        @inbounds active_node = nd[node_num]
        if (active_node.nodetype !== SUBEXPRESSION &&
            active_node.nodetype !== MOIVARIABLE &&
            active_node.nodetype !== VARIABLE &&
            active_node.nodetype !== VALUE)
            @inbounds children_idx = nzrange(adj, node_num)
            if (length(children_idx) > 0) # has any children
                for child in children_idx
                    is_match, match_num, match_dict = find_match(children_arr[child], nd, adj, const_list, parameter_values)
                    if ~is_match
                        @inbounds idx = children_arr[child]
                        @inbounds cn = nd[idx]
                        @inbounds parent_num = parent_dict[node_num]
                        push!(queue, (idx, node_num))
                        push!(new_nds, NodeData(cn.nodetype, cn.index, parent_num))
                        node_count += 1
                        parent_dict[idx] = node_count
                    else
                        node_count, const_len = substitute!(match_num, node_num, prior_prt,
                                                            nd, const_list, const_len, node_count,
                                                            parent_dict, match_dict, queue, new_nds)
                    end
                end
            end

        end
        println("(flatten iteration) end depth = $depth check")
    end
    expr.nd = new_nds
end
