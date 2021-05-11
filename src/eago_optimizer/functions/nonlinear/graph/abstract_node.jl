
"""
    NodeType

Each node in the directed graph can be classified into the following types
- VARIABLE: Denotes a decision variable.
- PARAMETER: An adjustable parameter value (not a decision variable).
- CONSTANT: A constant value
- EXPRESSION: Any other expression that isn't a subexpression
- SUBEXPRESSION: Any expression referencing a different graph representation.
"""
@enum(NodeClass, CONSTANT, PARAMETER, VARIABLE, EXPRESSION, SUBEXPRESSION)
#=
- SELECT: These nodes do not store values rather they reference the output
          stored in the buffer assocatied with an EXPRESSION. These nodes
          only occur after a multiple output function. In the mimo JuMP
          extension they correspond to intermediate variables introduced
          by the user.
=#

struct Variable end
struct Parameter end
struct Constant end
struct Expression end
struct Subexpression end
struct User end

function user end
function usern end

const ATOM_TYPE_INSTANCES = instances(AtomType)

abstract type AbstractNode end

"""
Describes connectivity and expression represented by node.
"""
struct Node <: AbstractNode
    node_class::NodeClass
    expr_type::AtomType
    first_index::Int
    second_index::Int
    arity::Int
    children::Vector{Int}
end

for (t, s, a) in ((Variable, VARIABLE, VAR_ATOM),
                  (Parameter, PARAMETER, PARAM_ATOM),
                  (Constant, CONSTANT, CONST_ATOM),)
    @eval function Node(::$t, i)
        return Node($s, $a, i, 0, 0, Int[])
    end
end

for (i, k) in enumerate(ALL_ATOM_TYPES)
    v = ALL_ATOM_DICT[k]
    @eval function Node(::typeof($v), children::Vector{Int})
        arity = length(children)
        return Node(EXPRESSION, $i, 0, arity, children)
    end
end

@inline _node_class(n::Node)   = n.node_class
@inline _expr_type(n::Node)    = n.node_expr_type
@inline _first_index(n::Node)  = n.first_index
@inline _second_index(n::Node) = n.node_second_index
@inline _arity(n::Node)        = n.arity
@inline _children(n::Node)     = n.children
@inline _child(n::Node, i)     = @inbounds getindex(n.children, i)

_create_call_node(nt, i, p, adj) = Node(ALL_ATOM_DICT[i], children)
function Node(d::JuMP._Derivatives.NodeData, c::AbstractVector{Int})
    nt = d.nodetype
    i = d.index
    p = d.parent
    if (nt == CALL) || (nt == CALLUNIVAR)
        return _create_call_node(nt, i, p, c)
    elseif nt == MOIVARIABLE
        error("MOI variable not supported.") # TODO: Confirm this error doesn't hit and delete.
    elseif nt == VALUE
        return Node(Value(), i)
    elseif nt == PARAMETER
        return Node(Parameter(), i)
    elseif nt == SUBEXPRESSION
        returnNode(Subexpression(), i)
    elseif nt == LOGIC
        error("Unable to load JuMP expression. Logical operators not currently supported.")
    elseif nt == COMPARISON
        error("Unable to load JuMP expression. Comparisons not currently supported.")
    end
    error("Node type not expected from JuMP.")
end

function _convert_node_list(x::Vector{JuMP._Derivatives.NodeData})
    y = Vector{Int}(undef, len)
    adj = adjmat(x)
    for i in eachindex(d)
        y[i] = Node(x[i], nzrange(adj, i))
    end
    return y
end
