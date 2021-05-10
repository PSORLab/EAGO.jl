
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
