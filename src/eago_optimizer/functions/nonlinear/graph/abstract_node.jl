
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
    ex_type::AtomType
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
    eval(quote
        function Node(::Val{true}, ::Val{$v}, c::Vector{Int})
            return Node(EXPRESSION, $v, 0, 0, length(c), c)
        end
    end)
end
@eval function Node(::Val{true}, ::Val{USERN}, i::Int, c::Vector{Int})
    return Node(EXPRESSION, USERN, i, 0, length(c), c)
end
for d in ALL_ATOM_TYPES
    @eval function Node(::Val{false}, ::Val{$d}, c::Vector{Int})
        return Node(EXPRESSION, $d, 0, 0, 1, c)
    end
end

_node_class(n::Node)   = n.node_class
_ex_type(n::Node)      = n.ex_type
_first_index(n::Node)  = n.first_index
_second_index(n::Node) = n.node_second_index
_arity(n::Node)        = n.arity
_children(n::Node)     = n.children
_child(n::Node, i)     = @inbounds getindex(n.children, i)

mv_eago_not_jump = setdiff(JuMP._Derivatives.operators,
                           union(Symbol[k for k in keys(REV_BIVARIATE_ATOM_DICT)],
                                 Symbol[k for k in keys(REV_NARITY_ATOM_DICT)]))
eago_mv_switch = quote end
for s in mv_eago_not_jump
    global eago_mv_switch = quote
        $eago_mv_switch
        (d == $s) && (return Node(Val(true), Val($s), v[c]))
    end
end
@eval function _create_call_node(i, v, c::UnitRange{Int}, op::OperatorRegistry)
    if i == 1
        return Node(Val(true), Val(PLUS), v[c])
    elseif i == 2
        return Node(Val(true), Val(MINUS), v[c])
    elseif i == 3
        return Node(Val(true), Val(MULT), v[c])
    elseif i == 4
        return Node(Val(true), Val(POW), v[c])
    elseif i == 5
        return Node(Val(true), Val(DIV), v[c])
    elseif i == 6
        error("If-else currently unsupported...")
    elseif i == 7
        return Node(Val(true), Val(MAX), v[c])
    elseif i == 8
        return Node(Val(true), Val(MIN), v[c])
    elseif i >= JuMP._Derivatives.USER_OPERATOR_ID_START
        i_mv = i - JuMP._Derivatives.USER_OPERATOR_ID_START + 1
        d = op.multivariate_id[i_mv]
        $eago_mv_switch
        return Node(Val(true), Val(USERN), i_mv, v[c])
    end
end

function binary_switch_typ(ids, exprs)
    if length(exprs) <= 3
        out = Expr(:if, Expr(:call, :(==), :i, ids[1]),
                   :(Node(Val(false), Val($(exprs[1])), v[c])))
        if length(exprs) > 1
            push!(out.args, binary_switch_typ(ids[2:end], exprs[2:end]))
        end
        return out
    else
        mid = length(exprs) >>> 1
        return Expr(:if, Expr(:call, :(<=), :i, ids[mid]),
                         binary_switch_typ(ids[1:mid], exprs[1:mid]),
                         binary_switch_typ(ids[mid+1:end], exprs[mid+1:end]))
    end
end

indx_JuMP = Int[]
indx_EAGO = AtomType[]
for k in univariate_operators
    if haskey(REV_UNIVARIATE_ATOM_DICT, k)
        k_EAGO = REV_UNIVARIATE_ATOM_DICT[k]
        push!(indx_JuMP, univariate_operator_to_id[k])
        push!(indx_EAGO, k_EAGO)
    end
end

uni_eago_not_jump = setdiff(univariate_operators, Symbol[k for k in keys(REV_UNIVARIATE_ATOM_DICT)])
uni_eago_not_jump = push!(uni_eago_not_jump, :-)
eago_uni_switch = quote end
for s in uni_eago_not_jump
    global eago_uni_switch = quote
        $eago_uni_switch
        (d == $s) && (return Node(Val(false), Val($s), v[c]))
    end
end
atom_switch = binary_switch_typ(indx_JuMP, indx_EAGO)
@eval function _create_call_node_uni(i::Int, v, c::UnitRange{Int}, op::OperatorRegistry)

    if i >= JuMP._Derivatives.USER_UNIVAR_OPERATOR_ID_START
        d = op.d.univariate_operator_to_id[i - JuMP._Derivatives.USER_UNIVAR_OPERATOR_ID_START + 1]
        $eago_uni_switch
        return Node(Val(true), Val(USER), i, v[c])
    end
    $atom_switch
end

function Node(d::JuMP._Derivatives.NodeData, child_vec, c::UnitRange{Int}, op::OperatorRegistry)
    nt = d.nodetype
    i = d.index
    if nt == JuMP._Derivatives.CALL
        return _create_call_node(i, child_vec, c, op)
    elseif nt == JuMP._Derivatives.CALLUNIVAR
        return _create_call_node_uni(i, child_vec, c, op)
    elseif nt == JuMP._Derivatives.VARIABLE
        return Node(Variable(), i)
    elseif nt == JuMP._Derivatives.MOIVARIABLE
        error("MOI variable not supported.") # TODO: Confirm this error doesn't hit and delete.
    elseif nt == JuMP._Derivatives.VALUE
        return Node(Constant(), i)
    elseif nt == JuMP._Derivatives.PARAMETER
        return Node(Parameter(), i)
    elseif nt == JuMP._Derivatives.SUBEXPRESSION
        return Node(Subexpression(), i)
    elseif nt == JuMP._Derivatives.LOGIC
        error("Unable to load JuMP expression. Logical operators not currently supported.")
    elseif nt == JuMP._Derivatives.COMPARISON
        error("Unable to load JuMP expression. Comparisons not currently supported.")
    end
    error("Node type = $nt not expected from JuMP.")
end

function _convert_node_list(x::Vector{JuMP._Derivatives.NodeData}, op::OperatorRegistry)
    y = Vector{Node}(undef, length(x))
    adj = JuMP._Derivatives.adjmat(x)
    child_vec = rowvals(adj)
    for i in eachindex(x)
        y[i] = Node(x[i], child_vec, nzrange(adj, i), op)
    end
    return y
end
