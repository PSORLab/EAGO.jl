
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
    @eval function Node(::$t, i::Int)
        return Node($s, $a, i, 0, 0, Int[])
    end
end

for v in (PLUS, MINUS, MULT, POW, DIV, MAX, MIN)
    @eval function Node(::Val{true}, ::Val{$v}, children::UnitRange{Int})
        arity = length(children)
        return Node(EXPRESSION, $v, 0, arity, children)
    end
end
@eval function Node(::Val{true}, ::Val{USERN}, i::Int, children::UnitRange{Int})
    arity = length(children)
    return Node(EXPRESSION, USERN, i, arity, children)
end
for d in ALL_ATOM_TYPES
    @eval function Node(::Val{false}, ::Val{$d}, children::Vector{Int})
        arity = length(children)
        return Node(EXPRESSION, $d, 0, 1, children)
    end
end

@inline _node_class(n::Node)   = n.node_class
@inline _expr_type(n::Node)    = n.node_expr_type
@inline _first_index(n::Node)  = n.first_index
@inline _second_index(n::Node) = n.node_second_index
@inline _arity(n::Node)        = n.arity
@inline _children(n::Node)     = n.children
@inline _child(n::Node, i)     = @inbounds getindex(n.children, i)

function _create_call_node(i, c::UnitRange{Int})
    if i == 1
        return Node(Val(true), Val(PLUS), c)
    elseif i == 2
        return Node(Val(true), Val(MINUS), c)
    elseif i == 3
        return Node(Val(true), Val(MULT), c)
    elseif i == 4
        return Node(Val(true), Val(POW), c)
    elseif i == 5
        return Node(Val(true), Val(DIV), c)
    elseif i == 6
        error("If-else currently unsupported...")
    elseif i == 7
        return Node(Val(Val(true), MAX), c)
    elseif i == 8
        return Node(Val(Val(true), MIN), c)
    elseif i >= JuMP._Derivatives.USER_OPERATOR_ID_START
        i_mv = i - JuMP._Derivatives.USER_OPERATOR_ID_START + 1
        return Node(Val(true), Val(USERN), i_mv, c)
    end
end
function _create_call_node_uni(i, c::UnitRange{Int})
    # TODO: Add binary switch here... for univariates
end

function Node(d::JuMP._Derivatives.NodeData, c::UnitRange{Int})
    nt = d.nodetype
    i = d.index
    if (nt == JuMP._Derivatives.CALL)
        return _create_call_node(i, c)
    elseif (nt == JuMP._Derivatives.CALLUNIVAR)
        return _create_call_node_uni(i, c)
    elseif nt == JuMP._Derivatives.MOIVARIABLE
        return Node(Variable(), i)
    elseif nt == JuMP._Derivatives.MOIVARIABLE
        error("MOI variable not supported.") # TODO: Confirm this error doesn't hit and delete.
    elseif nt == JuMP._Derivatives.VALUE
        return Node(Value(), i)
    elseif nt == JuMP._Derivatives.PARAMETER
        return Node(Parameter(), i)
    elseif nt == JuMP._Derivatives.SUBEXPRESSION
        return Node(Subexpression(), i)
    elseif nt == JuMP._Derivatives.LOGIC
        error("Unable to load JuMP expression. Logical operators not currently supported.")
    elseif nt == JuMP._Derivatives.COMPARISON
        error("Unable to load JuMP expression. Comparisons not currently supported.")
    end
    error("Node type not expected from JuMP.")
end

function _convert_node_list(x::Vector{JuMP._Derivatives.NodeData})
    y = Vector{Int}(undef, length(x))
    adj = JuMP._Derivatives.adjmat(x)
    for i in eachindex(x)
        y[i] = Node(x[i], nzrange(adj, i))
    end
    return y
end
