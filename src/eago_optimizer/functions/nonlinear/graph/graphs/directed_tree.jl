"""
    DirectedTree{S<:Real}

A tree graph with a single sink node.
"""
Base.@kwdef mutable struct DirectedTree{S<:Real} <: AbstractDirectedAcyclicGraph{S}
    "List of nodes"
    nodes::Vector{Node}                       = Node[]
    "List of index of variables"
    variables::Vector{Int}                    = Int[]
    "List of variable types"
    variable_types::Vector{VariableType}      = VariableType[]
    "List of constant values"
    constant_values::Vector{S}                = S[]
    "Number of nodes"
    node_count::Int                           = 0
    "Number of variables"
    variable_count::Int                       = 0
    "Number of constants"
    constant_count::Int                       = 0
    ""
    sparsity::Vector{Int}                     = Int[]
    ""
    rev_sparsity::Vector{Int}                 = Int[]
    dependent_variable_count::Int             = 0
    dependent_subexpression_count::Int        = 0
    dependent_subexpressions::Vector{Int}     = Int[]
    linearity::Linearity                      = CONSTANT
    user_operators::Ref{JuMP._Derivatives.UserOperatorRegistry}
end

const DAT = DirectedTree

# graph property access functions
@inbounds _nodes(g::DAT)           = g.nodes
@inbounds _variables(g::DAT)       = g.variables
@inbounds _variable_types(g::DAT)  = g.variable_types
@inbounds _constant_values(g::DAT) = g.constant_values

# user-define function access
@inline function _user_univariate_operator(g::DAT, i)
    @inbounds g.user_operators.univariate_operator_f[i]
end
@inline function _user_multivariate_operator(g::DAT, i)
    @inbounds g.user_operators.multivariate_operator_evaluator[i]
end


function fprop!(::Type{T}, g::DAT, b::AbstractCache) where {T<:AbstractCacheAttribute}
    f_init!(T, g, b)
    for k = length(g):-1:1
        nt = _node_type(g, k)
        if nt == EXPRESSION
            fprop!(T, Expression, g, b, k)
        elseif nt == VARIABLE
            fprop!(T, Variable, g, b, k)
        elseif nt == SUBEXPRESSION
            fprop!(T, Subexpression, g, b, k)
        end
    end
    return
end

function rprop!(::Type{T}, g::DAT, b::AbstractCache) where {T<:AbstractCacheAttribute}
    r_init!(T, g, b)
    for k = 1:length(g)
        nt = _node_type(g, k)
        if nt == EXPRESSION
            rprop!(T, Expression, g, b, k)
        elseif nt == VARIABLE
            rprop!(T, Variable, g, b, k)
        end
    end
    return
end

function DirectedTree{S}(sub::JuMP._SubexpressionStorage,
                         sub_sparsity::Dict{Int,Vector{Int}},
                         subexpr_linearity) where S<:Real

    nd = copy(func.nd)
    adj = copy(func.adj)
    const_values = copy(func.const_values)

    sparsity = copy(func.sparsity)
    for nd in func.nd
        if nd.nodetype === JuMP._Derivatives.SUBEXPRESSION
            append!(sparsity, sub_sparsity[nd.index])
        end
    end
    unique!(sparsity)
    sort!(sparsity)

    rev_sparsity_len = sparsity[end]
    rev_sparsity = zeros(Int, rev_sparsity_len)
    current_sparsity = grad_sparsity[1]
    current_sparsity_cnt = 1
    for i = 1:rev_sparsity_len
        if i == current_sparsity
            rev_sparsity[i] = current_sparsity_cnt
            current_sparsity_cnt += 1
            if current_sparsity_cnt <= length(grad_sparsity)
                current_sparsity = grad_sparsity[current_sparsity_cnt]
            else
                break
            end
        end
    end
    DirectedTree{S}(dependent_variable_count = length(grad_sparsity),
                    dependent_subexpressions = copy(func.dependent_subexpressions),
                    dependent_subexpression_count = length(dependent_subexpressions),
                    linearity = linearity(nd, adj, subexpr_linearity)
                    )
end
