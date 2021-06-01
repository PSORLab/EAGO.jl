
# Define standard forward and reverse propagation to switch of expression definitions
# for expressions.
function binary_switch(ids; is_forward = true)
    if length(ids) <= 3
        if is_forward
            out = Expr(:if, Expr(:call, :(==), :id, ids[1]),
                       :(return fprop!(t, Val($(ids[1])), g, c, k)))
        else
            out = Expr(:if, Expr(:call, :(==), :id, ids[1]),
                       :(return rprop!(t, Val$((ids[1])), g, c, k)))
        end
        if length(ids) > 1
            push!(out.args, binary_switch(ids[2:end]))
        end
        return out
    else
        mid = length(ids) >>> 1
        return Expr(:if, Expr(:call, :(<=), :id, ids[mid]),
                         binary_switch(ids[1:mid]),
                         binary_switch(ids[mid+1:end]))
    end
end

# Access gradient sparsity of JuMP storage.
_sparsity(d::JuMP._FunctionStorage) = d.grad_sparsity
_sparsity(d::JuMP._SubexpressionStorage) = d.sparsity

# Compute gradient sparsity from JuMP storage.
function _compute_sparsity(d::JuMP._FunctionStorage, sparse_dist::Dict{Int,Vector{Int}})
    sparsity = copy(_sparsity(d))
    for nd in d.nd
        if nd.nodetype === JuMP._Derivatives.SUBEXPRESSION
            append!(sparsity, sparse_dist[nd.index])
        end
    end
    unique!(sparsity)
    sort!(sparsity)
    sparsity, Int[]
end
function _compute_sparsity(d::JuMP._SubexpressionStorage, sparse_dist::Dict{Int,Vector{Int}})
    dep_subexpression = Int[]
    variable_dict = Dict{Int,Bool}()
    for n in d.nd
        if node.nodetype == JuMP._Derivatives.VARIABLE
            if !haskey(variable_dict, n.index)
                variable_dict[n.index] = true
            end
        end
        if node.nodetype == JuMP._Derivatives.SUBEXPRESSION
            push!(dep_subexpression, n.index)
        end
    end
    sparsity = collect(keys(variable_dict))
    unique!(dep_subexpressionession)
    sort!(dep_subexpression)
    sparsity, dep_subexpression
end

function linearity(d::JuMP._Derivatives.Linearity)
    if d == JuMP._Derivatives.LINEAR
        return LIN_LINEAR
    elseif d == JuMP._Derivatives.PIECEWISE_LINEAR
        return LIN_PIECEWISE_LINEAR
    elseif d == JuMP._Derivatives.NONLINEAR
        return LIN_NONLINEAR
    end
    # assumes d is then JuMP._Derivatives.CONSTANT
    return LIN_CONSTANT
end
function linearity(nd::Vector{JuMP._Derivatives.NodeData},
                   adj::SparseMatrixCSC{Bool,Int},
                   d::Vector{JuMP._Derivatives.Linearity})
    x = JuMP._Derivatives.classify_linearity(nd, adj, d)
    linearity.(x)
end

function OperatorRegistry(d::JuMP._Derivatives.UserOperatorRegistry)
    OperatorRegistry(Symbol[k for k in keys(d.multivariate_operator_to_id)],
                     d.multivariate_operator_to_id,
                     d.multivariate_operator_evaluator,
                     Symbol[k for k in keys(d.univariate_operator_to_id)],
                     d.univariate_operator_to_id,
                     d.univariate_operator_f,
                     d.univariate_operator_fprime,
                     d.univariate_operator_fprimeprime)
end
