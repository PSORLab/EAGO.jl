"""
    DirectedTree{S<:Real}

A tree graph with a single sink node.
"""
Base.@kwdef mutable struct DirectedTree{S<:Real} <: AbstractDirectedAcyclicGraph{S}
    "List of nodes"
    nodes::Vector{Node}                       = Node[]
    "List of index of variables in this tree"
    variables::Dict{Int,Int}                = Dict{Int,Int}()
    "Information on all variables..."
    v::VariableValues{S}                      = VariableValues{S}()
    "List of constant values"
    constant_values::Vector{S}                = S[]
    "Number of nodes"
    node_count::Int                           = 0
    "Number of variables"
    variable_count::Int                       = 0
    "Number of constants"
    constant_count::Int                       = 0
    sink_bnd::Interval{S}                     = Interval{S}(-Inf,Inf)
    ""
    sparsity::Vector{Int}                     = Int[]
    ""
    rev_sparsity::Dict{Int,Int}               = Dict{Int,Int}()
    dependent_variable_count::Int             = 0
    dep_subexpr_count::Int                    = 0
    dependent_subexpressions::Vector{Int}     = Int[]
    linearity::Linearity                      = LIN_CONSTANT
    user_operators::OperatorRegistry          = OperatorRegistry()
end
const DAT = DirectedTree

# graph property access functions
@inline _nodes(g::DAT)             = g.nodes
@inline _variables(g::DAT)         = g.variables
@inline _variable_types(g::DAT)    = g.v.variable_types
@inline _constant_values(g::DAT)   = g.constant_values
@inline _dep_subexpr_count(g::DAT) = g.dep_subexpr_count
# Each tree has a unique sparsity for a DAT since there is a single sink
@inline _sparsity(g::DAT, i)                  = g.sparsity
@inline _rev_sparsity(g::DAT, i::Int, k::Int) =  g.rev_sparsity[i]

# user-define function access
@inline function _user_univariate_operator(g::DAT, i)
    @inbounds g.user_operators.univariate_operator_f[i]
end
@inline function _user_multivariate_operator(g::DAT, i)
    @inbounds g.user_operators.multivariate_operator_evaluator[i]
end

function fprop!(t::T, g::DAT, b::AbstractCache) where {T<:AbstractCacheAttribute}
    f_init!(t, g, b)
    for k = _node_count(g):-1:1
        if _is_unlocked(b, k)
            nt = _node_class(g, k)
            if nt == EXPRESSION
                fprop!(t, Expression(), g, b, k)
            elseif nt == VARIABLE
                fprop!(t, Variable(), g, b, k)
            elseif nt == SUBEXPRESSION
                fprop!(t, Subexpression(), g, b, k)
            end
        end
    end
    return
end

function rprop!(t::T, g::DAT, b::AbstractCache) where {T<:AbstractCacheAttribute}
    flag = r_init!(t, g, b)
    for k = 1:_node_count(g)
        if _is_unlocked(b, k)
            nt = _node_class(g, k)
            if nt == EXPRESSION
                flag = rprop!(t, Expression(), g, b, k)
            elseif nt == VARIABLE
                flag = rprop!(t, Variable(), g, b, k)
            end
        end
    end
    return flag
end

# TODO Fix constructor...
function DirectedTree{S}(d, op::OperatorRegistry, sub_sparsity::Dict{Int,Vector{Int}}, subexpr_linearity) where S<:Real

    nd = copy(d.nd)
    adj = copy(d.adj)
    const_values = copy(d.const_values)

    sparsity, dependent_subexpressions = _compute_sparsity(d, sub_sparsity)
    rev_sparsity = Dict{Int,Int}()
    for (i,s) in enumerate(sparsity)
        rev_sparsity[s] = i
    end

    nodes = _convert_node_list(d.nd, op)
    lin = linearity(nd, adj, subexpr_linearity)
    DirectedTree{S}(nodes = nodes,
                    variables = rev_sparsity,
                    constant_values = const_values,
                    node_count = length(nodes),
                    variable_count = length(sparsity),
                    constant_count = length(const_values),
                    sparsity = sparsity,
                    rev_sparsity = rev_sparsity,
                    dependent_variable_count = length(sparsity),
                    dep_subexpr_count = length(dependent_subexpressions),
                    dependent_subexpressions = copy(dependent_subexpressions),
                    linearity = lin[1],
                    user_operators = op
                    )
end

#=
nd = wp._objective_nl.ex.nd
pushfirst!(nd, NodeData(JuMP._Derivatives.CALLUNIVAR, 2, -1))
nd[2] = NodeData(nd[2].nodetype, nd[2].index, 1)
for i = 3:length(nd)
    @inbounds nd[i] = NodeData(nd[i].nodetype, nd[i].index, nd[i].parent + 1)
end
I, J, V = findnz(wp._objective_nl.ex.adj)
I .+= 1
J .+= 1
pushfirst!(I, 2)
pushfirst!(J, 1)
pushfirst!(V, true)
m._working_problem._objective_nl.ex.adj = sparse(I, J, V)
=#

function _negate!(d::DirectedTree{S}) where S<:Real
    push!(d.nodes, Node(Val(false), Val(MINUS), Int[_node_count(d)]))
    d.node_count += 1
    d.sink_bnd = -d.sink_bnd
    return nothing
end
