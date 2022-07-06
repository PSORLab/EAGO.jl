
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

function Node(aux_info, d::JuMP._Derivatives.NodeData, child_vec, c::UnitRange{Int}, op)
    nt = d.nodetype
    i = d.index
    (nt == JuMP._Derivatives.CALL)          && return _create_call_node(i, child_vec, c, op)
    (nt == JuMP._Derivatives.CALLUNIVAR)    && return _create_call_node_uni(i, child_vec, c, op)
    (nt == JuMP._Derivatives.VALUE)         && return Node(Constant(), i)
    (nt == JuMP._Derivatives.PARAMETER)     && return Node(Parameter(), i)
    (nt == JuMP._Derivatives.SUBEXPRESSION) && return Node(Subexpression(), i)
    (nt == JuMP._Derivatives.VARIABLE)      && return !is_auxiliary_variable(aux_info, i) ? Node(Variable(), i) : Node(Select(), i)
    (nt == JuMP._Derivatives.LOGIC)         && error("Unable to load JuMP expression. Logical operators not currently supported.")
    (nt == JuMP._Derivatives.COMPARISON)    && error("Unable to load JuMP expression. Comparisons not currently supported.")
    error("Node type = $nt not expected from JuMP.")
end

function _convert_node_list(aux_info, x::Vector{JuMP._Derivatives.NodeData}, op)
    y = Vector{Node}(undef, length(x))
    adj = JuMP._Derivatives.adjmat(x)
    child_vec = rowvals(adj)
    for i in eachindex(x)
        y[i] = Node(aux_info, x[i], child_vec, nzrange(adj, i), op)
    end
    return y
end

# Access gradient sparsity of JuMP storage.
sparsity(d::JuMP._FunctionStorage) = d.grad_sparsity
sparsity(d::JuMP._SubexpressionStorage) = d.sparsity

# Compute gradient sparsity from JuMP storage.
function _compute_sparsity(d::JuMP._FunctionStorage, sparse_dict::Dict{Int,Vector{Int}}, is_sub, subexpr_indx)
    dep_subexpression = Int[]
    variable_dict = Dict{Int,Bool}()
    for n in d.nd
        if n.nodetype == JuMP._Derivatives.VARIABLE
            if !haskey(variable_dict, n.index)
                variable_dict[n.index] = true
            end
        end
        if n.nodetype == JuMP._Derivatives.SUBEXPRESSION
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
function _compute_sparsity(d::JuMP._SubexpressionStorage, sparse_dict::Dict{Int,Vector{Int}}, is_sub, subexpr_indx)
    dep_subexpression = Int[]
    variable_dict = Dict{Int,Bool}()
    for n in d.nd
        if n.nodetype == JuMP._Derivatives.VARIABLE
            if !haskey(variable_dict, n.index)
                variable_dict[n.index] = true
            end
        end
        if n.nodetype == JuMP._Derivatives.SUBEXPRESSION
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

function linearity(d::JuMP._Derivatives.Linearity)
    (d == JuMP._Derivatives.LINEAR)            && return LIN_LINEAR
    (d == JuMP._Derivatives.PIECEWISE_LINEAR)  && return LIN_PIECEWISE_LINEAR
    (d == JuMP._Derivatives.NONLINEAR)         && return LIN_NONLINEAR
    LIN_CONSTANT                               # assumes d is then JuMP._Derivatives.CONSTANT
end
function linearity(nd::Vector{JuMP._Derivatives.NodeData}, adj::SparseMatrixCSC{Bool,Int}, d::Vector{JuMP._Derivatives.Linearity})
    x = JuMP._Derivatives.classify_linearity(nd, adj, d)
    linearity.(x)
end

function OperatorRegistry(d::JuMP._Derivatives.UserOperatorRegistry)
    mv_id = collect(keys(d.multivariate_operator_to_id))
    mv_operator_to_id = d.multivariate_operator_to_id
    mv_operator_evaluator = d.multivariate_operator_evaluator
    u_operator_id = collect(keys(d.univariate_operator_to_id))
    u_to_id = d.univariate_operator_to_id
    u_to_f = d.univariate_operator_f
    u_fprime = d.univariate_operator_fprime
    u_fprimeprime = d.univariate_operator_fprimeprime
    OperatorRegistry(mv_id, mv_operator_to_id, mv_operator_evaluator, u_operator_id, u_to_id, u_to_f, u_fprime, u_fprimeprime)
end