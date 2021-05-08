
"""
    NodeType

Each node in the directed graph can be classified into the following types
- VARIABLE: Denotes a decision variable.
- PARAMETER: An adjustable parameter value (not a decision variable).
- CONSTANT: A constant value
- EXPRESSION: Any other expression that isn't a subexpression
- SUBEXPRESSION: Any expression referencing a different graph representation.
"""
@enum(NodeType, CONSTANT, PARAMETER, VARIABLE, EXPRESSION, SUBEXPRESSION)
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
    node_type::NodeType
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

for (k,v) in ALL_ATOM_DICT
    @eval function Node(::typeof($v), children::Vector{Int})
        arity = length(children)
        return Node(EXPRESSION, $k, 0, arity, children)
    end
end

@inline _node_type(n::Node)    = n.node_type
@inline _expr_type(n::Node)    = n.node_expr_type
@inline _first_index(n::Node)  = n.first_index
@inline _second_index(n::Node) = n.node_second_index
@inline _arity(n::Node)        = n.arity
@inline _children(n::Node)     = n.children
@inline _child(n::Node, i)     = @inbounds getindex(n.children, i)
