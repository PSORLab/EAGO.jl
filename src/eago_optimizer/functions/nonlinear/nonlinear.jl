# Copyright (c) 2018: Matthew Wilhelm & Matthew Stuber.
# This work is licensed under the Creative Commons Attribution-NonCommercial-
# ShareAlike 4.0 International License. To view a copy of this license, visit
# http://creativecommons.org/licenses/by-nc-sa/4.0/ or send a letter to Creative
# Commons, PO Box 1866, Mountain View, CA 94042, USA.
#############################################################################
# EAGO
# A development environment for robust and global optimization
# See https://github.com/PSORLab/EAGO.jl
#############################################################################
# TODO
#############################################################################


"""
$(FUNCTIONAME)

Stores a general quadratic function with a buffer.
"""
mutable struct NonlinearExpression{V} <: AbstractEAGOConstraint

    "List of nodes in nonlinear expression"
    nd::Vector{JuMP.NodeData}
    "Adjacency Matrix for the expression"
    adj::SparseMatrixCSC{Bool,Int64}
    const_values::Vector{Float64}
    affine_term::SAF

    setstorage::Vector{V}
    numberstorage::Vector{Float64}
    isnumber::Vector{Bool}
    value::V
    value_available::Bool

    tp1storage::Vector{Float64}
    tp2storage::Vector{Float64}
    tp3storage::Vector{Float64}
    tp4storage::Vector{Float64}
    tpdict::Dict{Int64,Tuple{Int64,Int64,Int64,Int64}}

    # sparsity of constraint + indices in node to reference
    grad_sparsity::Vector{Int64}       # indices of variables in the problem space (size = np)
    is_branch::Vector{Bool}            # is variable (size = np1)
    branch_indices::Vector{Int64}      # map variable pos to pos in node (-1 if nonbranching)

    # role in problem
    dependent_variable_count::Int
    dependent_subexpression_count::Int
    dependent_subexpressions::Vector{Int64}
    linearity::JuMP._Derivatives.Linearity
end

"""
$(FUNCTIONAME)

Stores a general nonlinear function with a buffer represented by the sum of a tape
and a scalar affine function.
"""
mutable struct BufferedNonlinearFunction{V} <: AbstractEAGOConstraint
    expr::NonlinearExpression{V}
    saf::SAF
    lower_bound::Float64
    upper_bound::Float64
    last_relax_convex::Bool
    last_relax_concave::Bool
    last_past_reverse::Bool
end

###
### Constructor utilities
###

function extract_affine_term!(d::NonlinearExpression{V}) where V
    # TODO
    return nothing
end

function extract_affine_term!(d::BufferedNonlinearFunction{V}) where V
    extract_affine_term!(d.expr)
    return nothing
end

###
### Constructor definitions
###
function NonlinearExpression{V}(sub::JuMP._SubexpressionStorage,
                                     subexpr_linearity::Vector{JuMP._Derivatives.Linearity},
                                     tag::T) where T
    nd = copy(func.nd)
    adj = copy(func.adj)
    const_values = copy(func.const_values)

    setstorage = fill(MC{N,T}(Interval(-Inf, Inf)), lenx)
    numberstorage = zeros(lenx)
    isnumber = fill(false, lenx)

    tpdict = Dict{Int64,Tuple{Int64,Int64,Int64,Int64}}()
    tp1_count = 0
    tp2_count = 0
    for i = 1:lenx
        node = @inbounds nd[i]
        op = node.index
        if double_tp(op)
            tp1_count += 1
            tpdict[i] = (tp1_count, tp1_count, tp2_count, tp2_count)
        elseif single_tp(op)
            tp1_count += 1
            tp2_count += 1
            tpdict[i] = (tp1_count, tp1_count, -1, -1)
        end
    end
    tp1storage = zeros(tp1_count)
    tp2storage = zeros(tp1_count)
    tp3storage = zeros(tp2_count)
    tp4storage = zeros(tp2_count)

    dependent_subexpressions = copy(func.dependent_subexpressions)
    dependent_subexpression_count = length(dependent_subexpressions)

    linearity = JuMP._Derivatives.classify_linearity(nd, adj, subexpr_linearity)

    grad_sparsity = copy(func.grad_sparsity)  # sorted by JUmp, _FunctionStorage
    is_branch = Bool[]                        # set in label_branch_variables routine
    branch_indices = Bool[]                   # set in label_branch_variables routine

    grad_sparsity = Float64[]
    subexpression = NonlinearExpression{MC{N,T}}(nd, adj, const_values, SAF(SAT[], 0.0),
                                                            setstorage, numberstorge, isnumber, zero(MC{N,T}),
                                                            dependent_variable_count, dependent_subexpressions,
                                                            tp1storage, tp2storage, tp3storage, tp4storage, tpdict,
                                                            grad_sparsity, is_branch, branch_indices,
                                                            JuMP._Derivatives.CONSTANT)
    extract_affine_term!(subexpression)
end

function NonlinearExpression{V}()
    return BufferedNonlinearSubexpression{V}(JuMP.NodeData[], spzeros(Bool, 1), Float64[], SAF(SAT[], 0.0),
                                             V[], Float64[], Bool[], zero(V),
                                             Float64[], Float64[], Float64[], Float64[],
                                             Dict{Int64,Tuple{Int64,Int64,Int64,Int64}}()
                                             Int64[], Bool[], Int64[], 0, Int64[], JuMP._Derivatives.CONSTANT)
end

function BufferedNonlinearFunction(func::JuMP._FunctionStorage, bnds::MOI.NLPBoundsPair,
                           subexpr_linearity::Vector{JuMP._Derivatives.Linearity},
                           tag::T) where T <: RelaxTag

    nd = copy(func.nd)
    adj = copy(func.adj)

    const_values = copy(func.const_values)

    setstorage = fill(MC{N,T}(Interval(-Inf, Inf)), lenx)
    numberstorage = zeros(lenx)
    isnumber = fill(false, lenx)

    tpdict = Dict{Int64,Tuple{Int64,Int64,Int64,Int64}}()
    tp1_count = 0
    tp2_count = 0
    for i = 1:lenx
        node = @inbounds nd[i]
        op = node.index
        if double_tp(op)
            tp1_count += 1
            tpdict[i] = (tp1_count, tp1_count, tp2_count, tp2_count)
        elseif single_tp(op)
            tp1_count += 1
            tp2_count += 1
            tpdict[i] = (tp1_count, tp1_count, -1, -1)
        end
    end
    tp1storage = zeros(tp1_count)
    tp2storage = zeros(tp1_count)
    tp3storage = zeros(tp2_count)
    tp4storage = zeros(tp2_count)

    dependent_variable_count = length(func.grad_sparsity)
    saf_buffer = SAF(SAT[SAT(0.0, VI(-1)) for i=1:dependent_variable_count], 0.0)

    # sorted by JUmp, _FunctionStorage
    grad_sparsity = copy(func.grad_sparsity)

    # set in label_branch_variables routine
    is_branch = Bool[]
    branch_indices = Bool[]

    dependent_subexpressions = copy(func.dependent_subexpressions)
    dependent_subexpression_count = length(dependent_subexpressions)

    linearity = JuMP._Derivatives.classify_linearity(nd, adj, subexpr_linearity)

    nonlinear_constraint =  BufferedNonlinearFunction{MC{N,T}}(nd, adj, const_values,  SAF(SAT[], 0.0),
                                                       setstorage, numberstorage, isnumber, saf_buffer, zero(MC{N,T}),
                                                       tp1storage, tp2storage, tp3storage, tp4storage, tpdict,
                                                       grad_sparsity, is_branch, branch_indices,
                                                       dependent_variable_count, dependent_subexpression_count,
                                                       dependent_subexpressions, bnds, linearity)
    return extract_affine_term(nonlinear_constraint)
end

function BufferedNonlinearFunction{V}()
    return BufferedNonlinearFunction{V}(NonlinearExpression{V}(), SAF(SAT[], 0.0), MOI.NLPBoundsPair(-Inf, Inf))
end

function set_node_flag!(f::BufferedNonlinearFunction{V}) where V
    d.has_value = false
    d.last_relax_convex = false
    d.last_relax_concave = false

    return nothing
end

"""
$(TYPEDEF)

Extracts the `convex` affine relaxaiton is `use_cvx` to `f.saf` then adds the `affine_terms` to
this to form the affine relaxation of the function.
"""
function unpack_value!(f::BufferedNonlinearFunction{V}, use_cvx::Bool) where V

    return nothing
end

###
### Defines evaluator storage structure
###
"""
$(TYPEDEF)

MOI.AbstractNLPEvaluator for calculating relaxations of nonlinear terms.

Checks that the resulting value should be a number...

$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct Evaluator <: MOI.AbstractNLPEvaluator

    user_operators::JuMP._Derivatives.UserOperatorRegistry = JuMP._Derivatives.UserOperatorRegistry()
    has_user_mv_operator::Bool = false
    parameter_values::Vector{Float64} = Float64[]

    current_node::NodeBB = NodeBB()
    lower_variable_bounds::Vector{Float64} = Float64[]
    upper_variable_bounds::Vector{Float64} = Float64[]
    x_value::Vector{Float64} = Float64[]
    ni_map::Vector{Int64}

    "Context used to guard against domain violations & branch on these violations if necessary"
    subgrad_tighten::Bool = false
    subgrad_tighten_reverse::Bool = false
    ctx::GuardCtx = GuardCtx()

    subexpressions::Vector{NonlinearExpression}
    subexpressions_eval::Vector{Bool}
end

"""
$(FUNCTIONNAME)

Sets the current node in the Evaluator structure.
"""
function set_node!(evaluator::Evaluator, n::NodeBB)
    evaluator.current_node = NodeBB(n)
    for i = 1:length(evaluator.current_node)
        @inbounds evaluator.lower_variable_bounds[ni_map[i]] = n.lower_variable_bounds[i]
        @inbounds evaluator.upper_variable_bounds[ni_map[i]] = n.lower_variable_bounds[i]
    end
    fill!(evaluator.subexpressions_eval, false)

    return nothing
end

function set_reference_point!(evaluator::Evaluator, x::Vector{Float64})
    fill!(evaluator.subexpressions_eval, false)

    return nothing
end

function retrieve_node(d::Evaluator)
    cn = d.current_node
    return NodeBB(copy(d.lower_variable_bounds[ni_map]),
                  copy(d.upper_variable_bounds[ni_map]),
                  cn.lower_bound, cn.upper_bound, cn.depth, cn.id)
end

# Returns false if subexpression has been evaluated at current reference point
prior_eval(d::Evaluator, i::Int64) = @inbounds subexpressions_eval[i]

#=
Assumes the sparsities are sorted...
=#
function copy_subexpression_value!(k::Int, op::Int, setstorage::Vector{MC{N1,T}}, substorage::MC{N2,T},
                                   cv_buffer::Vector{Float64}, cc_buffer::Vector{Float64},
                                   func_sparsity::Vector{Int64}, sub_sparsity::Vector{Int64}) where {N1, N2, T <: RelaxTag}

    sset = @inbounds #TODO

    # fill cv_grad/cc_grad buffers
    fill!(cv_buffer, 0.0)
    fill!(cc_buffer, 0.0)

    sub_sparsity_count = 1
    subs_index = @inbounds sub_sparsity[1]
    for i = 1:length(func_sparsity)
        func_index = @inbounds func_sparsity[i]
        if func_index === subs_index
            @inbounds cv_buffer[i] = sset.cv_grad[sub_sparsity_count]
            @inbounds cc_buffer[i] = sset.cc_grad[sub_sparsity_count]
            sub_sparsity_count += 1
            subs_index = @inbounds sub_sparsity[sub_sparsity_count]
        end
    end

    cv_grad = SVector(cv_buffer)
    cc_grad = SVector(cc_buffer)

    setstorage[k] = MC{N1,T}(sset.cv, sset.cc, sset.Intv, cv_grad, cc_grad, sset.cnst)

    return nothing
end

include("forward_pass.jl")

###
### Define forward evaluation pass
###
function forward_pass!(evaluator::Evaluator, d::NonlinearExpression{V}) where V
    # check that prior subexpressions have been evaluated
    # i.e. box_id is same and reference point is the same
    for i = 1:d.dependent_subexpression_count
        if !prior_eval(evaluator, i)
            forward_pass!(evaluator, subexpr)
        end
    end
    forward_pass_kernel!(d.nd, d.adj, d.x, evaluator.lower_variable_bounds,
                         evaluator.upper_variable_bounds, d.setstorage,
                         d.numberstorage, d.numvalued, d.tpdict,
                         d.tp1storage, d.tp2storage, d.tp3storage, d.tp4storage,
                         evaluator.user_operators, evaluator.subexpression_isnum,
                         evaluator.subexpr_values_flt, evaluator.subexpr_values_set,
                         evaluator.num_mv_buffer, d.set_mv_buffer, evaluator.ctx,
                         evaluator.is_post, evaluator.is_intersect,
                         evaluator.is_first_eval, d.cv_grad_buffer, d.cc_grad_buffer)
    return nothing
end

function forward_pass!(evaluator::Evaluator, d::BufferedNonlinearFunction{V}) where V
    forward_pass!(evaluator, d.expr)
    set_value!(d.expr, d.expr.value ∩ d.bnds)
    d.has_value = true
    d.last_past_reverse = false
    return nothing
end

###
### Define backwards evaluation pass
###

include("reverse_pass.jl")

"""
$(FUNCTIONNAME)

A reverse_pass! on a `BufferedNonlinear` structure `d` intersects the existing value of the `d` with
constraint bounds then reverse propagates a set-valued operator (by default McCormick operator) along the
computational tape. The tapes are updated in place and boolean value is returned indicating whether the
reverse propagation yeilded a infeasible point (true = still feasible, false is proved infeasible).
"""
function reverse_pass!(evaluator::Evaluator, d::NonlinearExpression{V}) where V
    return reverse_pass_kernel!(d.nd, d.adj, d.x, evaluator.lower_variable_bounds,
                                evaluator.upper_variable_bounds, d.setstorage,
                                d.numberstorage, d.numvalued, evaluator.subexpression_isnum,
                                evaluator.subexpr_values_set, evaluator.is_post)
end

function reverse_pass!(evaluator::Evaluator, d::BufferedNonlinearFunction{V}) where V
    d.last_past_reverse = true
    set_value!(d.expr, d.expr.value ∩ d.bnds)
    return reverse_pass!(evaluator, d.expr)
end

###
### Interval bounding definitions
###
function lower_interval_bound(d::BufferedNonlinearFunction{V}, n::NodeBB) where V
    if !d.has_value
        forward_pass!(d.evaluator, d)
    end
    return get_lo(get_value(d))
end

function interval_bound(d::BufferedNonlinearFunction{V}, n::NodeBB) where V
    if !d.has_value
        forward_pass!(d.evaluator, d)
    end
    return get_interval(get_value(d))
end

###
### Parsing definitions
###
function eliminate_fixed_variables!(f::NonlinearExpression{V}}, v::Vector{VariableInfo}) where V
    num_constants = length(f.const_values)
    indx_to_const_loc = Dict{Int,Int}()
    for i = 1:length(expr.nd)
        nd = @inbounds expr.nd[i]
        # Assumes MOI Variable have been eliminated previously...
        if nd.nodetype === JuMP._Derivatives.VARIABLE
            indx = nd.index
            if v[indx].is_fixed
                if haskey(indx_to_const_loc, indx)
                    const_loc = indx_to_const_loc[indx]
                    expr.nd[i] = NodeData(JuMP._Derivatives.VALUE, const_loc, nd.parent)
                else
                    push!(const_values, v[indx].lower_bound)
                    num_constants += 1
                    indx_to_const_loc[indx] = num_constants
                    expr.nd[i] = NodeData(nd.nodetype, num_constants, nd.parent)
                end
                f.isnumber[i] = true
            end
        end
    end

    return nothing
end

function eliminate_fixed_variables!(f::BufferedNonlinearFunction{V}, v::Vector{VariableInfo}) where V
    eliminate_fixed_variables!(f.expr, v)
end
