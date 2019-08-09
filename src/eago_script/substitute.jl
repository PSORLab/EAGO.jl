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
end
function Template_Graph(nd::Dict{Int,Template_Node}, dag::Vector{Pair{Int,Int}})
    Template_Graph([nd[i] for i in 1:length(nd)], dag, length(nd), length(dag))
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

function matching_info(x::Template_Node, y::NodeData)
    match_flag = false
    if (x.type == :op)
        if (y.nodetype == CALL)
            if (operator_to_id[x.value] == y.index)
                match_flag = true
            end
        elseif (y.nodetype == CALLUNIVAR)
            if univariate_operator_to_id[x.value] == y.index
                match_flag = true
        #elseif y.nodetype == SUBEXPRESSION TODO: Add subexpression handling latter
        end
    elseif (x.type == :num)
        if x.check(x.num_value)
            match_flag = true
        end
    elseif (x.type == :expr)
        match_flag = true
    else
        error("invalid type symbol")
    end
end

function is_match(pattern::Template_Graph, indx::Int, nd::Vector{NodeData}, dag_adj::SparseMatrixCSC{Bool,Int})

    match_flag = true

    # make pattern adjacency matrix
    pattern_length = pattern.ndlen
    dag_length = pattern.daglen
    adj = spzeros(Bool, pattern_length, pattern_length)
    for i in 1:dag_length
        @inbounds p = pattern.dag[i]
        @inbounds x = p[1]
        @inbounds y = p[2]
        adj[x, y] = true
    end
    pat_children_arr = rowvals(adj)
    dag_children_arr = rowvals(dag_adj)

    # do a breadth first search with paired template, nd data,
    # if any pair of children fail then
    #symbol_lib = Dict{:Symbol, }
    pindx_initial = 1
    queue = Pair{Int,Int}[(pindx_initial, indx)]
    while ~isempty(queue) && (match_flag == true)
        (num_pat, num_dag) = popfirst!(queue)
        if matching_info(pattern.nd[num_pat], nd[indx])
            @inbounds pat_children_idx = nzrange(adj, num_pat)
            pat_length = length(pat_children_idx)
            if pat_length > 0
                @inbounds dag_children_idx = nzrange(dag_adj, num_dag)
                if pat_length == length(dag_children_idx)
                    # at this point children should match and
                    # if symbol is registered... then expressions should be equal
                else
                    match_flag = false
                    break
                end
            else
                match_flag = false
                break
            end
        else
            match_flag = false
            break
        end
    end

    return match_flag
end

function find_match(indx::Int, nd::Vector{NodeData}, adj::SparseMatrixCSC{Bool,Int})
    flag = false
    pattern_number = -1
    @inbounds sub_len = DAG_LENGTHS[1]
    for i in 1:sub_len
        @inbounds pattern = DAG_PATTERNS[i]
        if is_match(pattern, indx, nd, adj)
            flag = true
            pattern_number = i
            break
        end
    end
    return flag, pattern_number
end

# we assume a tree structure, so if we don't load child nodes,
# then they are effectively deleted
function substitute!(match_num::Int, node_num::Int, nd::Vector{NodeData},
                     parent_dict::Dict{Int,Int}, queue::Vector{Tuple{Int,Int}},
                     new_nds::Vector{NodeData})

    match_template = DAG_PATTERNS[match_num]
    subst_template = DAG_SUBSTITUTIONS[match_num]
    # create correspondence between pattern and substitutions
    # add nodes to graph... adding any terminal children nodes to queue
end

# searchs through expression breadth first search that short-cirucits
# if match found
#
function flatten_expression!(expr::_NonlinearExprData)
    nd = expr.nd
    adj = adjmat(nd)
    children_arr = rowvals(adj)
    node_count = 1

    # if the first node matches a pattern then substitute... else load the first node
    @inbounds active_node = nd[1]
    if (active_node.nodetype !== SUBEXPRESSION &&
        active_node.nodetype !== MOIVARIABLE &&
        active_node.nodetype !== VARIABLE &&
        active_node.nodetype !== VALUE)
        is_match, match_num = find_match(1, nd, adj)
        if is_match
            node_count += substitute!(match_num, node_num, nd, parent_dict, queue, new_nds)
        else
            new_nds = NodeData[active_node]
            parent_dict[1] = 1
            push!(queue, (1,-1))
        end
    end

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
                    is_match, match_num = findmatch(children_arr[child], nd, adj)
                    if ~is_match
                        @inbounds idx = children_arr[child]
                        @inbounds cn = nd[idx]
                        @inbounds parent_num = parent_dict[node_num]
                        push!(queue, (idx, node_num))
                        push!(new_nds, NodeData(cn.nodetype, cn.index, parent_num))
                        node_count += 1
                        parent_dict[idx] = node_count
                    else
                        node_count += substitute!(match_num, node_num, nd, parent_dict, queue, new_nds)
                    end
                end
            end

        end
    end
    expr.nd = new_nds
end
