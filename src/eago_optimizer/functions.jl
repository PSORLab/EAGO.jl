"""
$(TYPEDEF)

An abstract super-type used for representing constraints built by EAGO's backend.
"""
abstract type AbstractEAGOConstraint end

"""

Computes the lower interval bound for `AbstractEAGOConstraint` representing an
inequality constraint.
"""
function lower_interval_bound end

"""

Computes a tuple representing the lower and upper interval bounds for a
`AbstractEAGOConstraint` representing an equality constraint.
"""
function interval_bound end


####
#### Affine Storage
####

"""
$(TYPEDEF)

"""
struct AffineFunctionIneq <: AbstractEAGOConstraint
    terms::Vector{Tuple{Float64,Int}}
    constant::Float64
    len::Int
end

function lower_interval_bound(f::AffineFunctionIneq, y::NodeBB)

    terms = f.terms
    lo_bnds = y.lower_variable_bounds
    up_bnds = y.upper_variable_bounds

    lower_interval_bound = f.constant
    for i = 1:f.len
        coeff, indx = @inbounds terms[k]
        if coeff > 0.0
            lvb = @inbounds lo_bnds[indx]
            lower_interval_bound += coeff*lvb
        else
            uvb = @inbounds up_bnds[indx]
            lower_interval_bound += coeff*uvb
        end
    end

    return lower_interval_bound
end


"""
$(TYPEDEF)

"""
struct AffineFunctionEq <: AbstractEAGOConstraint
    terms::Vector{Tuple{Float64,Int}}
    constant::Float64
end

function interval_bound(f::AffineFunctionEq, y::NodeBB)
    terms = f.terms
    lo_bnds = y.lower_variable_bounds
    up_bnds = y.upper_variable_bounds

    lower_interval_bound = f.constant
    upper_interval_bound = f.constant
    for i = 1:f.len
        coeff, indx = @inbounds terms[k]
        lvb = @inbounds lo_bnds[indx]
        uvb = @inbounds up_bnds[indx]
        if coeff > 0.0
            lower_interval_bound += coeff*lvb
            upper_interval_bound += coeff*uvb
        else
            lower_interval_bound += coeff*uvb
            upper_interval_bound += coeff*lvb
        end
    end

    return lower_interval_bound, upper_interval_bound
end


####
#### Quadratic Storage
####
"""
$(TYPEDEF)

Stores a general quadratic function with a buffer.
"""
mutable struct BufferedQuadraticIneq <: AbstractEAGOConstraint
    func::SQF
    buffer::OrderedDict{Int, Float64}
    saf::SAF
    nx::Int
end

function lower_interval_bound(f::BufferedQuadraticIneq, n::NodeBB)

    lo_bnds = n.lower_variable_bounds
    up_bnds = n.upper_variable_bounds
    lower_interval_bound = Interval{Float64}(f.func.constant)

    for aff_term in f.func.affine_terms
        coeff = aff_term.coefficient
        vi = aff_term.variable_index.value
        @inbounds xL = lo_bnds[vi]
        @inbounds xU = up_bnds[vi]
        lower_interval_bound += coeff > 0.0 ? coeff*xL : coeff*xU
    end

    for quad_term in f.func.quadratic_terms
        coeff = quad_term.coefficient
        vi1 = quad_term.variable_index_1.value
        vi2 = quad_term.variable_index_2.value
        @inbounds xL = lo_bnds[vi1]
        @inbounds xU = up_bnds[vi1]
        if vi1 === vi2
            if coeff > 0.0
                lower_interval_bound += (0.0 < xL) ? coeff*xL*xL : ((xU <= 0.0) ? coeff*xU*xU : 0.0)
            else
                lower_interval_bound += (xL < xU) ? coeff*xU*xU : coeff*xL*xL
            end
        else
            @inbounds il2b = lo_bnds[vi2]
            @inbounds iu2b = up_bnds[vi2]
            lower_interval_bound += coeff*Interval{Float64}(xL, xU)*Interval{Float64}(il2b, iu2b)
        end
    end

    return val_intv.lo, val_intv.hi
end

"""
$(TYPDEF)

Stores a general quadratic function with a buffer.
"""
mutable struct BufferedQuadraticEq <: AbstractEAGOConstraint
    func::SQF
    buffer::OrderedDict{Int, Float64}
    saf::SAF
    nx::Int
end

function interval_bound(f::BufferedQuadraticEq, n::NodeBB)

    lo_bnds = n.lower_variable_bounds
    up_bnds = n.upper_variable_bounds
    val_intv = Interval(f.func.constant)

    for aff_term in f.func.affine_terms
        coeff = aff_term.coefficient
        vi = aff_term.variable_index.value
        @inbounds il1b = lo_bnds[vi]
        @inbounds iu1b = up_bnds[vi]
        val_intv += coeff*Interval(il1b, iu1b)
    end

    for quad_term in f.func.quadratic_terms
        coeff = quad_term.coefficient
        vi1 = quad_term.variable_index_1.value
        vi2 = quad_term.variable_index_2.value
        @inbounds il1b = lo_bnds[vi1]
        @inbounds iu1b = up_bnds[vi1]
        if vi1 === vi2
            val_intv += coeff*pow(Interval(il1b, iu1b), 2)
        else
            @inbounds il2b = lo_bnds[vi2]
            @inbounds iu2b = up_bnds[vi2]
            val_intv += coeff*Interval(il1b, iu1b)*Interval(il2b, iu2b)
        end
    end

    return val_intv.lo, val_intv.hi
end

#=
NodeType convention is defined to parallel JuMP's nonlinear evaluator
@enum NodeType: CALL CALLUNIVAR VARIABLE VALUE SUBEXPRESSION PARAMETER

const OPERATORS = [:+, :-, :*, :^, :/, :ifelse, :max, :min]
const USER_OPERATOR_ID_START = length(operators) + 1
const OPERATOR_TO_ID = Dict{Symbol,Int}()
for i = 1:length(OPERATORS)
    OPERATOR_TO_ID[OPERATORS[i]] = i
end

const UNIVARIATE_OPERATORS = Symbol[:+, :-, :abs]
const UNIVARIATE_OPERATOR_TO_ID = Dict{Symbol,Int}(:+ => 1, :- => 2, :abs => 3)
const UNIVARIATE_OPERATOR_DERIV = Any[:(one(x)), :(-one(x)), :(ifelse(x >= 0, one(x), -one(x)))]
for (op, deriv) in Calculus.symbolic_derivatives_1arg()
    push!(UNIVARIATE_OPERATORS, op)
    push!(UNIVARIATE_OPERATOR_DERIV, deriv)
    UNIVARIATE_OPERATOR_TO_ID[op] = length(UNIVARIATE_OPERATORS)
end
const USER_UNIVAR_OPERATOR_ID_START = length(UNIVARIATE_OPERATORS) + 1

"""
$(TYPDEF)

Stores a general nonlinear function with a buffer.
"""
mutable struct BufferedNonlinearEq{N, T<:RelaxTag} <: AbstractEAGOConstraint
    "Main evaluator"
    evaluator::Evaluator
    "List of nodes in nonlinear expression"
    node_list::NodeData
    const_values::Vector{Float64}
    set_storage::Vector{MC{N,T}}
    grad_sparsity::Vector{Int64}
    dependent_subexpressions::Vector{Int64}
end

function interval_bound(f::BufferedNonlinearEq{N,T}, n::NodeBB) where {N, T<:RelaxTag}
end

function relax!(m::Optimizer, f::BufferedNonlinearEq{N,T}, indx::Int, check_safe::Bool) where {N, T<:RelaxTag}
end

"""
$(FUNCTIONAME)

Stores a general quadratic function with a buffer.
"""
mutable struct BufferedNonlinear <: AbstractEAGOConstraint
    func
end
=#

####
#### Nonlinear Storage
####

struct EmptyNLPEvaluator <: MOI.AbstractNLPEvaluator
    _current_node::NodeBB
    has_nlobj::Bool
end
EmptyNLPEvaluator() = EmptyNLPEvaluator(NodeBB(),false)
set_current_node!(x::EmptyNLPEvaluator, n::NodeBB) = ()

MOI.features_available(::EmptyNLPEvaluator) = [:Grad, :Jac, :Hess]
MOI.initialize(::EmptyNLPEvaluator, features) = nothing
MOI.eval_objective(::EmptyNLPEvaluator, x) = NaN
function MOI.eval_constraint(::EmptyNLPEvaluator, g, x)
    @assert length(g) == 0
    return
end
MOI.eval_objective_gradient(::EmptyNLPEvaluator, g, x) = nothing
MOI.jacobian_structure(::EmptyNLPEvaluator) = Tuple{Int64,Int64}[]
MOI.hessian_lagrangian_structure(::EmptyNLPEvaluator) = Tuple{Int64,Int64}[]
function MOI.eval_constraint_jacobian(::EmptyNLPEvaluator, J, x)
    @assert length(J) == 0
    return
end
function MOI.eval_hessian_lagrangian(::EmptyNLPEvaluator, H, x, σ, μ)
    @assert length(H) == 0
    return
end

empty_nlp_data() = MOI.NLPBlockData([], EmptyNLPEvaluator(), false)
