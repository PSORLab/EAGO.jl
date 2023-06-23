
# Define standard forward and reverse propagation to switch of expression definitions
# for expressions.
function binary_switch(ids; is_forward = true)
    if length(ids) <= 3
        if is_forward
            out = Expr(:if, Expr(:call, :(==), :id, ids[1]), :(return fprop!(t, Val($(ids[1])), g, c, k)))
        else
            out = Expr(:if, Expr(:call, :(==), :id, ids[1]), :(return rprop!(t, Val$((ids[1])), g, c, k)))
        end
        (length(ids) > 1) && push!(out.args, binary_switch(ids[2:end]))
        return out
    else
        mid = length(ids) >>> 1
        return Expr(:if, Expr(:call, :(<=), :id, ids[mid]), binary_switch(ids[1:mid]), binary_switch(ids[mid+1:end]))
    end
end

function Node(aux_info, d::MOINL.Node, child_vec, c::UnitRange{Int}, op)
    nt = d.type
    i = d.index
    (nt == MOINL.NODE_CALL_MULTIVARIATE)    && return _create_call_node(i, child_vec, c, op)
    (nt == MOINL.NODE_CALL_UNIVARIATE)      && return _create_call_node_uni(i, child_vec, c, op)
    (nt == MOINL.NODE_VALUE)                && return Node(Constant(), i)
    (nt == MOINL.NODE_PARAMETER)            && return Node(Parameter(), i)
    (nt == MOINL.NODE_SUBEXPRESSION)        && return Node(Subexpression(), i)
    (nt == MOINL.NODE_VARIABLE)             && return !is_auxiliary_variable(aux_info, i) ? Node(Variable(), i) : Node(Select(), i)
    (nt == MOINL.NODE_MOI_VARIABLE)         && return !is_auxiliary_variable(aux_info, i) ? Node(Variable(), i) : Node(Select(), i)
    (nt == MOINL.NODE_LOGIC)                && error("Unable to load JuMP expression. Logical operators not currently supported.")
    (nt == MOINL.NODE_COMPARISON)           && error("Unable to load JuMP expression. Comparisons not currently supported.")
    error("Node type = $nt not expected from JuMP.")
end

function _convert_node_list(aux_info, x::Vector{MOINL.Node}, op)
    y = Vector{Node}(undef, length(x))
    adj = MOINL.adjacency_matrix(x)
    child_vec = rowvals(adj)
    for i in eachindex(x)
        y[i] = Node(aux_info, x[i], child_vec, nzrange(adj, i), op)
    end
    return y
end

# Access gradient sparsity of MOI storage.
sparsity(d::MOIRAD._FunctionStorage) = d.grad_sparsity

# Compute gradient sparsity from MOI storage.
function _compute_sparsity(d::MOIRAD._FunctionStorage, sparse_dict::Dict{Int,Vector{Int}}, is_sub, subexpr_indx)
    dep_subexpression = Int[]
    variable_dict = Dict{Int,Bool}()
    for n in d.nodes
        if n.type == MOINL.NODE_VARIABLE
            if !haskey(variable_dict, n.index)
                variable_dict[n.index] = true
            end
        end
        if n.type == MOINL.NODE_SUBEXPRESSION
            push!(dep_subexpression, n.index)
        end
    end
    sparsity = collect(keys(variable_dict))
    unique!(dep_subexpression)
    sort!(dep_subexpression)
    for s in dep_subexpression
        append!(sparsity, sparse_dict[s])
    end
    unique!(sparsity)
    sort!(sparsity)
    if is_sub
        sparse_dict[subexpr_indx] = sparsity
    end
    sparsity, dep_subexpression
end
function _compute_sparsity(d::MOIRAD._SubexpressionStorage, sparse_dict::Dict{Int,Vector{Int}}, is_sub, subexpr_indx)
    dep_subexpression = Int[]
    variable_dict = Dict{Int,Bool}()
    for n in d.nodes
        if n.type == MOINL.NODE_VARIABLE
            if !haskey(variable_dict, n.index)
                variable_dict[n.index] = true
            end
        end
        if n.type == MOINL.NODE_SUBEXPRESSION
            push!(dep_subexpression, n.index)
        end
    end
    sparsity = collect(keys(variable_dict))
    unique!(dep_subexpression)
    sort!(dep_subexpression)
    for s in dep_subexpression
        append!(sparsity, sparse_dict[s])
    end
    unique!(sparsity)
    sort!(sparsity)
    if is_sub
        sparse_dict[subexpr_indx] = sparsity
    end
    sparsity, dep_subexpression
end

function linearity(d::MOIRAD.Linearity)
    (d == MOIRAD.LINEAR)            && return LIN_LINEAR
    (d == MOIRAD.PIECEWISE_LINEAR)  && return LIN_PIECEWISE_LINEAR
    (d == MOIRAD.NONLINEAR)         && return LIN_NONLINEAR
    LIN_CONSTANT                   # assumes d is then MOINL.CONSTANT
end

function linearity(nd::Vector{MOINL.Node}, adj::SparseMatrixCSC{Bool,Int}, d::Vector{MOIRAD.Linearity})
    x = MOIRAD._classify_linearity(nd, adj, d)
    linearity.(x)
end

function OperatorRegistry(d::MOINL.OperatorRegistry)
    u_operators = d.univariate_operators
    u_operator_id = collect(keys(d.univariate_operator_to_id))
    u_to_id = d.univariate_operator_to_id
    u_operator_start = d.univariate_user_operator_start
    r_u_operators = d.registered_univariate_operators
    mv_operators = d.multivariate_operators
    mv_operator_id = collect(keys(d.multivariate_operator_to_id))
    mv_to_id = d.multivariate_operator_to_id
    mv_operator_start = d.multivariate_user_operator_start
    r_mv_operators = d.registered_multivariate_operators
    l_operators = d.logic_operators
    l_operator_id = collect(keys(d.logic_operator_to_id))
    l_to_id = d.logic_operator_to_id
    c_operators = d.comparison_operators
    c_operator_id = collect(keys(d.comparison_operator_to_id))
    c_to_id = d.comparison_operator_to_id
    OperatorRegistry(u_operators, u_operator_id, u_to_id, u_operator_start, r_u_operators, mv_operators, mv_operator_id, mv_to_id, mv_operator_start, r_mv_operators, l_operators, l_operator_id, l_to_id, c_operators, c_operator_id, c_to_id)
end