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

"""
$(TYPEDEF)

A structure used to hold objectives and constraints added to EAGO model.
The constraints generally aren't used for relaxations.

$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct InputProblem

    # variables
    _variable_info::Vector{VariableInfo} = VariableInfo[]
    _variable_count::Int64 = 0

    # last constraint index added
    _last_constraint_index::Int = 0

    # linear constraints
    _linear_leq_constraints::Vector{Tuple{SAF, LT}} = Tuple{SAF, LT}[]
    _linear_geq_constraints::Vector{Tuple{SAF, GT}} = Tuple{SAF, GT}[]
    _linear_eq_constraints::Vector{Tuple{SAF, ET}} = Tuple{SAF, ET}[]

    _linear_leq_count::Int = 0
    _linear_geq_count::Int = 0
    _linear_eq_count::Int = 0

    # quadratic constraints
    _quadratic_leq_constraints::Vector{Tuple{SQF, LT}} = Tuple{SQF, LT}[]
    _quadratic_geq_constraints::Vector{Tuple{SQF, GT}} = Tuple{SQF, GT}[]
    _quadratic_eq_constraints::Vector{Tuple{SQF, ET}} = Tuple{SQF, ET}[]

    _quadratic_leq_count::Int = 0
    _quadratic_geq_count::Int = 0
    _quadratic_eq_count::Int = 0

    # conic constraints
    _conic_second_order::Vector{Tuple{VECOFVAR, MOI.SecondOrderCone}} = Tuple{VECOFVAR, MOI.SecondOrderCone}[]

    _conic_second_order_count::Int = 0

    # objectives
    _objective_sv::SV = SV(-1)
    _objective_saf::SAF = SAF(SAT[], 0.0)
    _objective_sqf::SQF = SQF(SAT[], SQT[], 0.0)
    _objective_type::ObjectiveType = UNSET

    # nlp constraints
    _nlp_data::MOI.NLPBlockData = empty_nlp_data()

    # attributes
    _optimization_sense::MOI.OptimizationSense = MOI.MIN_SENSE
end

Base.@kwdef mutable struct VariableNodeMap
    eq_variable_indx::OrderDict{CI{SV, ET}, Int} = OrderDict{CI{SV, ET}, Int}()
    leq_variable_indx::OrderDict{CI{SV, LT}, Int} = OrderDict{CI{SV, LT}, Int}()
    geq_variable_indx::OrderDict{CI{SV, GT}, Int} = OrderDict{CI{SV, GT}, Int}()
end
getindex(v::VariableNodeMap, i::CI{SV, ET}) = v.eq_variable_indx[i]
getindex(v::VariableNodeMap, i::CI{SV, LT}) = v.leq_variable_indx[i]
getindex(v::VariableNodeMap, i::CI{SV, GT}) = v.geq_variable_indx[i]

setindex!(v::VariableNodeMap, i::CI{SV, ET}, val) = (v.eq_variable_indx[i] = val)
setindex!(v::VariableNodeMap, i::CI{SV, LT}, val) = (v.leq_variable_indx[i] = val)
setindex!(v::VariableNodeMap, i::CI{SV, GT}, val) = (v.geq_variable_indx[i] = val)

#=
Holds specialized constraint functions used by EAGO to generate cuts
=#
Base.@kwdef mutable struct ParsedProblem

    # Problem classification
    _problem_type::ProblemType = UNCLASSIFIED

    # objectives
    _objective_sv::SV = SV(-1)
    "_objective_saf stores the objective and is used for constructing linear affine cuts
     of any ObjectiveType"
    _objective_saf::SAF = SAF(SAT[], 0.0)
    _objective_sqf::BufferedQuadraticIneq = BufferedQuadraticIneq()
    _objective_nl = nothing
    _objective_type::ObjectiveType = UNSET

    # constraints
    _saf_leq::Vector{AffineFunctionIneq} = Vector{AffineFunctionIneq}[]
    _saf_eq::Vector{AffineFunctionEq} = Vector{AffineFunctionEq}[]
    _sqf_leq::Vector{BufferedQuadraticIneq} = Vector{BufferedQuadraticIneq}[]
    _sqf_eq::Vector{BufferedQuadraticEq} = Vector{BufferedQuadraticEq}[]

    _saf_leq_count::Int = 0
    _saf_eq_count::Int = 0
    _sqf_leq_count::Int = 0
    _sqf_eq_count::Int = 0

    # nlp constraints
    _nlp_data::MOI.NLPBlockData = empty_nlp_data()

    # variables
    _variable_info::Vector{VariableInfo} = VariableInfo[]
    _variable_count::Int = 0

    _var_leq_count::Int = 0
    _var_geq_count::Int = 0
    _var_eq_count::Int = 0

    # need to retreive primal _relaxed_variable_index
    _relaxed_variable_node_map::VariableNodeMap
    _relaxed_variable_index::Vector{VI} = VI[]
    _relaxed_variable_eq::Vector{CI{SV, ET}} = CI{SV, ET}[]
    _relaxed_variable_lt::Vector{CI{SV, LT}} = CI{SV, LT}[]
    _relaxed_variable_gt::Vector{CI{SV, GT}} = CI{SV, GT}[]
end

function bound_objective(m::Optimizer)
    d = x._relaxed_evaluator

    if x._objective_type === NONLINEAR

        objective_lo = eval_objective_lo(d)
        constraints = d.constraints
        constr_num = d.constraint_number
        constraints_intv_lo = zeros(Float64, constr_num)
        constraints_intv_hi = zeros(Float64, constr_num)
        eval_constraint_lo!(d, constraints_intv_lo)
        eval_constraint_hi!(d, constraints_intv_hi)
        constraints_bnd_lo = d.constraints_lbd
        constraints_bnd_hi = d.constraints_ubd

        for i = 1:d.constraint_number
            @inbounds constraints_intv_lo = constraints_bnd_lo[i]
            @inbounds constraints_intv_hi = constraints_bnd_hi[i]
            if (constraints_intv_lo > constraints_intv_hi) || (constraints_intv_hi < constraints_intv_lo)
                feas = false
                break
            end
        end
    elseif x._objective_type ===  SINGLE_VARIABLE
        obj_indx = x._objective_sv.variable.value
        objective_lo = @inbounds y.lower_variable_bounds[obj_indx]
    elseif x._objective_type === SCALAR_AFFINE
        objective_lo = interval_bound(x._objective_saf, y, true)
    elseif x._objective_type === SCALAR_QUADRATIC
        objective_lo = interval_bound(x._objective_sqf, y, true)
    end

    if objective_lo > x._lower_objective_value
        x._lower_objective_value = objective_lo
    end
end
