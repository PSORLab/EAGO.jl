
# Define standard forward and reverse propagation to switch of expression definitions
# for expressions.
function binary_switch(ids, exprs; is_rev = true)
    if length(exprs) <= 3
        if is_rev
            out = Expr(:if, Expr(:call, :(==), :id, ids[1]),
                       :(fprop!(T, $(exprs[1]), g, b, k)))
        else
            out = Expr(:if, Expr(:call, :(==), :id, ids[1]),
                       :(rprop!(T, $(exprs[1]), g, b, k)))
        end
        if length(exprs) > 1
            push!(out.args, binary_switch(ids[2:end], exprs[2:end]))
        end
        return out
    else
        mid = length(exprs) >>> 1
        return Expr(:if, Expr(:call, :(<=), :id, ids[mid]),
                         binary_switch(ids[1:mid], exprs[1:mid]),
                         binary_switch(ids[mid+1:end], exprs[mid+1:end]))
    end
end


# Access gradient sparsity of JuMP storage.
_sparsity(d::JuMP._FunctionStorage) where S<:Real = d.grad_sparsity
_sparsity(d::JuMP._SubexpressionStorage) where S<:Real = d.sparsity

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
    for (i,n) in enumerate(d.nd)
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
