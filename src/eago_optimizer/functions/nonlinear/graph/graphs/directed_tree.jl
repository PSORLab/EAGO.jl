"""
    AbstractCache

Abstract supertype used for information storage object the directed acyclic graph.
"""
abstract type AbstractCache end

"""
    AbstractCacheAttribute

Abstract supertype used for attributes stored in a cache.
"""
abstract type AbstractCacheAttribute end

Base.@kwdef mutable struct VariableValues{T<:Real}
    x0::Vector{T}                          = T[]
    x::Vector{T}                           = T[]
    lower_variable_bounds::Vector{T}       = T[]
    upper_variable_bounds::Vector{T}       = T[]
    node_to_variable_map::Vector{Int}      = Int[]
    variable_to_node_map::Vector{Int}      = Int[]
    variable_types::Vector{VariableType}   = VariableType[]
end

@inline _val(b::VariableValues{T}, i::Int) where T = @inbounds b.x[i]
@inline _lbd(b::VariableValues{T}, i::Int) where T = @inbounds b.lower_variable_bounds[i]
@inline _ubd(b::VariableValues{T}, i::Int) where T = @inbounds b.upper_variable_bounds[i]
_val(b::VariableValues{T}) where T = b.x
_lbd(b::VariableValues{T}) where T = b.lower_variable_bounds
_ubd(b::VariableValues{T}) where T = b.upper_variable_bounds

function _get_x!(::Type{BranchVar}, out::Vector{T}, v::VariableValues{T}) where T<:Real
    @inbounds for i = 1:length(v.node_to_variable_map)
        out[i] = v.x[v.node_to_variable_map[i]]
    end
    return nothing
end

function _initialize_or_copy!(y::VariableValues{T}, x::VariableValues{T}, s::Symbol) where {T<:Real}
    isempty(getfield(y, s)) ? setfield!(y, s, copy(getfield(x, s))) : copy!(getfield(y, s), getfield(x, s))
end
function update_box_and_pnt!(y::VariableValues{T}, x::VariableValues{T}, update_box::Bool) where {T<:Real}
    if update_box
        _initialize_or_copy!(y, x, :x0)
        _initialize_or_copy!(y, x, :lower_variable_bounds)
        _initialize_or_copy!(y, x, :upper_variable_bounds)
    end
    _initialize_or_copy!(y, x, :x)
end

"""
    DirectedTree

A tree graph with a single sink node.
"""
Base.@kwdef mutable struct DirectedTree <: AbstractDirectedAcyclicGraph
    "List of nodes"
    nodes::Vector{Node}                       = Node[]
    "List of index of variables in this tree"
    variables::Dict{Int,Int}                  = Dict{Int,Int}()
    "Information on all variables..."
    v::VariableValues{Float64}                = VariableValues{Float64}()
    "List of constant values"
    constant_values::Vector{Float64}          = Float64[]
    "List of constant values"
    parameter_values::Vector{Float64}         = Float64[]
    "Number of nodes"
    node_count::Int                           = 0
    "Number of variables"
    variable_count::Int                       = 0
    "Number of constants"
    constant_count::Int                       = 0
    sink_bnd::Interval{Float64}               = Interval{Float64}(-Inf,Inf)
    ""
    sparsity::Vector{Int}                     = Int[]
    ""
    rev_sparsity::Dict{Int,Int}               = Dict{Int,Int}()
    dependent_variable_count::Int             = 0
    dep_subexpr_count::Int                    = 0
    dependent_subexpressions::Vector{Int}     = Int[]
    linearity::Linearity                      = LIN_CONSTANT
    user_operators::OperatorRegistry          = OperatorRegistry()
    #children::SpraseMatrixCSC{Bool,Int}      = spzeros(Bool,1,1)
    #parents::SparseMatrixCSC{Bool,Int}       = spzeros(Bool,1,1)
end
#DirectedTree(n::Int) = DirectedTree(children = spzeros(Bool,n,n), parents=spzeros(Bool,n,n))
const DAT = DirectedTree

# node property access functions that can be defined at abstract type
_node(g::DAT, i)            = g.nodes[i]
_nodes(g::DAT)              = g.nodes

_variable(g::DAT, i)        = g.variables[i]
_variables(g::DAT)          = g.variables

_constant_value(g::DAT, i)  = g.constant_values[i]
_constant_values(g::DAT)    = g.constant_values
_parameter_value(g::DAT, i) = g.parameter_values[i]
_parameter_values(g::DAT)   = g.parameter_values

_node_class(g::DAT, i)      = _node_class(_node(g, i))
_ex_type(g::DAT, i)         = _ex_type(_node(g, i))
_first_index(g::DAT, i)     = _first_index(_node(g, i))
_secondary_index(g::DAT, i) = _secondary_index(_node(g, i))
_arity(g::DAT, i)           = _arity(_node(g, i))
_children(g::DAT, i)        = _children(_node(g, i))
_child(g::DAT, i, j)        = _child(_node(g, j), i)

_node_count(g::DAT)     = g.node_count
_variable_count(g::DAT) = g.variable_count
_constant_count(g::DAT) = g.constant_count

_dep_subexpr_count(g::DAT)            = length(g.dependent_subexpressions)
_sparsity(g::DAT, i)                  = g.sparsity
_rev_sparsity(g::DAT, i::Int, k::Int) = g.rev_sparsity[i]

_user_univariate_operator(g::DAT, i) = g.user_operators.univariate_operator_f[i]
_user_multivariate_operator(g::DAT, i) = g.user_operators.multivariate_operator_evaluator[i]


function DirectedTree(aux_info, d, op::OperatorRegistry, sub_sparsity::Dict{Int,Vector{Int}}, subexpr_linearity, parameter_values, is_sub, subexpr_indx)

    nd = copy(d.nd)
    adj = copy(d.adj)
    const_values = copy(d.const_values)

    sparsity, dependent_subexpressions = _compute_sparsity(d, sub_sparsity, is_sub, subexpr_indx)
    rev_sparsity = Dict{Int,Int}()
    for (i,s) in enumerate(sparsity)
        rev_sparsity[s] = i
    end

    nodes = _convert_node_list(aux_info, d.nd, op)
    lin = linearity(nd, adj, subexpr_linearity)
    DirectedTree(nodes = nodes,
                    variables = rev_sparsity,
                    constant_values = const_values,
                    parameter_values = parameter_values,
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

forward_uni = [i for i in instances(AtomType)]
setdiff!(forward_uni, [VAR_ATOM; PARAM_ATOM; CONST_ATOM; SELECT_ATOM; SUBEXPR])
f_switch = binary_switch(forward_uni, is_forward = true)
@eval function fprop!(t::T, ex::Expression, g::DAT, c::AbstractCache , k::Int) where T<:AbstractCacheAttribute
    id = _ex_type(g, k)
    $f_switch
    error("fprop! for ex_type = $id not defined.")
    return
end

reverse_uni = [i for i in instances(AtomType)]
setdiff!(reverse_uni, [VAR_ATOM; PARAM_ATOM; CONST_ATOM; SELECT_ATOM; SUBEXPR])
r_switch = binary_switch(reverse_uni, is_forward = false)
@eval function rprop!(t::T, ex::Expression, g::DAT, c::AbstractCache, k::Int) where T<:AbstractCacheAttribute
    id = _ex_type(g, k)
    $r_switch
    error("rprop! for ex_type = $id not defined.")
    return
end