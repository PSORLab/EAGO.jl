# Copyright (c) 2018: Matthew Wilhelm & Matthew Stuber.
# This code is licensed under MIT license (see LICENSE.md for full details)
#############################################################################
# EAGO
# A development environment for robust and global optimization
# See https://github.com/PSORLab/EAGO.jl
#############################################################################
# src/eago_optimizer/constraints.jl
# Defines constraints supported by optimizer and how to store them.
#############################################################################

# Sets used in general constraints
const INEQ_SETS = Union{LT, GT, ET}
const VAR_SETS = Union{LT, GT, ET, ZO, MOI.Integer}

##### Utilities for checking that JuMP model contains variables used in expression
function check_inbounds!(m::Optimizer, vi::VI)
    if !(1 <= vi.value <= m._input_problem._variable_count)
        error("Invalid variable index $vi. ($(m._input_problem._variable_count) variables in the model.)")
    end
    return nothing
end
check_inbounds!(m::Optimizer, f::SAF) = foreach(x -> check_inbounds!(m, x.variable), f.terms)
check_inbounds!(m::Optimizer, f::VECOFVAR) = foreach(x -> check_inbounds!(m, x), f.variables)
function check_inbounds!(m::Optimizer, f::SQF)
    foreach(x -> check_inbounds!(m, x.variable), f.affine_terms)
    for term in f.quadratic_terms
        check_inbounds!(m, term.variable_1)
        check_inbounds!(m, term.variable_2)
    end
    return nothing
end

MOI.supports_constraint(::Optimizer, ::Type{VI}, ::Type{S}) where S <: VAR_SETS = true
MOI.supports_constraint(::Optimizer,::Type{T},::Type{S}) where {T<:Union{SAF,SQF},S<:INEQ_SETS} = true

MOI.is_valid(m::Optimizer, v::VI) = (1 <= v.value <= m._input_problem._variable_count)
MOI.is_valid(m::Optimizer, c::CI{VI,S}) where S <: VAR_SETS = (1 <= c.value <= m._input_problem._variable_count)
MOI.is_valid(m::Optimizer, c::CI{F,S}) where {F<:Union{SAF,SQF}, S<:INEQ_SETS} = (1 <= c.value <= length(_constraints(m,F,S)))

MOI.get(m::Optimizer, ::MOI.NumberOfConstraints{VI,S}) where S<:VAR_SETS = length(_constraints(m,VI,S))
MOI.get(m::Optimizer, ::MOI.NumberOfConstraints{F,S}) where {F<:Union{SAF,SQF},S<:INEQ_SETS} = length(_constraints(m,F,S))

MOI.get(m::Optimizer, ::MOI.ListOfConstraintIndices{VI,S}) where S<:VAR_SETS = collect(keys(_constraints(m,VI,S)))
MOI.get(m::Optimizer, ::MOI.ListOfConstraintIndices{F,S}) where {F<:Union{SAF,SQF},S<:INEQ_SETS} = collect(keys(_constraints(m,F,S)))

MOI.add_variable(m::Optimizer) = VI(m._input_problem._variable_count += 1)

function MOI.add_constraint(m::Optimizer, f::F, s::S) where {F<:Union{SAF,SQF},S<:INEQ_SETS}
    check_inbounds!(m, f)
    ci = CI{F, S}(m._input_problem._constraint_count += 1)
    _constraints(m, F, S)[ci] = (f, s)
    return ci
end
function MOI.add_constraint(m::Optimizer, f::VI, s::S) where S<:VAR_SETS
    check_inbounds!(m, f)
    ci = CI{VI,S}(f.value)
    _constraints(m,VI,S)[ci] = (f, s)
    return ci
end

result_index_1_error(v::T) where T = throw(MOI.ResultIndexBoundsError{T}(v, 1))
function MOI.get(model::Optimizer, v::MOI.VariablePrimal, vi::MOI.VariableIndex)
    check_inbounds!(model, vi)
    (v.result_index != 1) && result_index_1_error(v)
    return model._global_optimizer._continuous_solution[vi.value]
end
function MOI.get(m::Optimizer{R,S,T}, v::MOI.ConstraintPrimal, c::CI{VI,<:Any}) where {R,S,T}
    (v.result_index != 1) && result_index_1_error(v)
    return MOI.get(m, MOI.VariablePrimal(), MOI.VariableIndex(c.value))
end
function MOI.get(m::Optimizer{R,Q,T}, v::MOI.ConstraintPrimal, c::CI{F,S}) where {R,Q,T,F,S}
    (v.result_index != 1) && result_index_1_error(v) 
    return _constraint_primal(m._global_optimizer, F, S)[c]
end

MOI.get(m::Optimizer, ::MOI.ConstraintFunction, c::CI{F,S}) where {F,S} = _constraints(m,F,S)[c][1]
MOI.get(m::Optimizer, ::MOI.ConstraintSet, c::CI{F,S}) where {F,S} = _constraints(m,F,S)[c][2]

function MOI.empty!(m::Optimizer{R,S,T}) where {R,S,T}

    MOI.empty!(m.subsolver_block)
    MOI.empty!(m._global_optimizer)
    m._input_problem = InputProblem()
    m._working_problem = ParsedProblem()

    m._termination_status_code = MOI.OPTIMIZE_NOT_CALLED
    m._result_status_code = MOI.OTHER_RESULT_STATUS
    m._run_time = 0.0
    m._objective_value  = -Inf
    m._objective_bound  =  Inf
    m. _relative_gap     = Inf
    m._iteration_count  = 0
    m._node_count       = 0

    return nothing
end

function MOI.is_empty(m::Optimizer{R,S,T}) where {R,S,T}
                           
    flag = true
    flag &= MOI.is_empty(m._global_optimizer)
    flag &= isempty(m._input_problem)
    flag &= isempty(m._working_problem)
    flag &= isempty(m.subsolver_block)
    flag &= m._termination_status_code == MOI.OPTIMIZE_NOT_CALLED
    flag &= m._result_status_code == MOI.OTHER_RESULT_STATUS

    # set constructor reset on empty! and to zero in initial_parse! in parse.jl
    flag &= iszero(m._run_time)
    flag &= iszero(m._iteration_count)
    flag &= iszero(m._node_count)
    flag &= m._objective_value == -Inf
    flag &= m._objective_bound ==  Inf
    flag &= m. _relative_gap   == Inf

    return flag
end

MOI.supports_incremental_interface(m::Optimizer) = true
MOI.copy_to(model::Optimizer, src::MOI.ModelLike) = MOIU.default_copy_to(model, src)

#####
##### Set & get attributes of model
#####
MOI.supports(::Optimizer, ::MOI.Silent) = true

function MOI.set(m::Optimizer, s::MOI.Silent, value)
    push!(m._optimizer_attributes_set, s)
    unique!(m._optimizer_attributes_set)
    if value
        m._parameters.verbosity = 0
        m._parameters.log_on = false
    else
        m._parameters.verbosity = 1
    end
    return
end

function MOI.set(m::Optimizer, s::MOI.TimeLimitSec, ::Nothing)
    push!(m._optimizer_attributes_set, s)
    unique!(m._optimizer_attributes_set)
    m._parameters.time_limit = Inf
end
function MOI.set(m::Optimizer, s::MOI.TimeLimitSec, v::Int)
    push!(m._optimizer_attributes_set, s)
    unique!(m._optimizer_attributes_set)
    m._parameters.time_limit = v
end
function MOI.set(m::Optimizer, s::MOI.TimeLimitSec, v::Float64)
    push!(m._optimizer_attributes_set, s)
    unique!(m._optimizer_attributes_set)
    m._parameters.time_limit = v
end

MOI.get(m::Optimizer, ::MOI.ListOfOptimizerAttributesSet) = m._optimizer_attributes_set

function MOI.get(m::Optimizer, ::MOI.ListOfConstraintTypesPresent)
    constraint_types = []
    for S in (ZO, MOI.Integer)
        if MOI.get(m, MOI.NumberOfConstraints{VI,S}()) > 0
            push!(constraint_types, (VI,S))
        end
    end
    for S in (LT, GT, ET), T in (VI, SAF, SQF)
        if MOI.get(m, MOI.NumberOfConstraints{T,S}()) > 0
            push!(constraint_types, (VI,S))
        end
    end
    return constraint_types
end

MOI.get(m::Optimizer, v::MOI.ObjectiveValue) = !isone(v.result_index) ? result_index_1_error(v) : m._objective_value
MOI.get(m::Optimizer, v::MOI.PrimalStatus) = !isone(v.result_index) ? MOI.NO_SOLUTION : m._result_status_code
MOI.get(m::Optimizer, ::MOI.DualStatus) = MOI.NO_SOLUTION
MOI.get(m::Optimizer, ::MOI.ObjectiveBound) = m._objective_bound
MOI.get(m::Optimizer, ::MOI.NumberOfVariables) = m._input_problem._variable_count
MOI.get(m::Optimizer, ::MOI.SolverName) = "EAGO: Easy Advanced Global Optimization"
MOI.get(m::Optimizer, ::MOI.SolverVersion) = "0.8.1"
MOI.get(m::Optimizer, ::MOI.TerminationStatus) = m._termination_status_code
MOI.get(m::Optimizer, ::MOI.SolveTimeSec) = m._run_time
MOI.get(m::Optimizer, ::MOI.NodeCount) = m._node_count
MOI.get(m::Optimizer, ::MOI.ResultCount) = (m._result_status_code === MOI.FEASIBLE_POINT) ? 1 : 0
MOI.get(m::Optimizer, ::MOI.TimeLimitSec) = m._parameters.time_limit
MOI.get(m::Optimizer, ::MOI.Silent) = m._parameters.verbosity == 0
MOI.get(m::Optimizer, ::MOI.ListOfVariableIndices) = [VI(i) for i = 1:m._input_problem._variable_count]

function MOI.get(m::Optimizer, ::MOI.RelativeGap)
    b = MOI.get(m, MOI.ObjectiveBound())
    v = MOI.get(m, MOI.ObjectiveValue())
    return relative_gap(b,v)
end

_to_sym(d) = error("EAGO only supports raw parameters with Symbol or String names.")
_to_sym(d::String) = Symbol(d)
_to_sym(d::Symbol) = d
function MOI.get(m::Optimizer, p::MOI.RawOptimizerAttribute)
    s = _to_sym(p.name)
    s in EAGO_PARAMETERS ? getfield(m._parameters, s) : getfield(m, s)
end

raw_param_name(p::MOI.RawOptimizerAttribute) = _to_sym(p.name)
raw_param_name(p) = nothing
function MOI.set(m::Optimizer, p::MOI.RawOptimizerAttribute, x)
    inds = findall(x -> raw_param_name(x) == raw_param_name(p), m._optimizer_attributes_set)
    deleteat!(m._optimizer_attributes_set, inds)
    push!(m._optimizer_attributes_set, p)
    s = _to_sym(p.name)
    if (s == :relaxed_optimizer) || (s == :upper_optimizer)
        setfield!(m, s, Incremental(x))
    else
        s in EAGO_PARAMETERS ? setfield!(m._parameters, s, x) : setfield!(m, s, x)
    end
end
MOI.get(m::Optimizer, p::MOI.RawStatusString) = string(m._global_optimizer._end_state)

#####
##### Support, set, and evaluate objective functions
#####
MOI.supports(::Optimizer, ::MOI.TimeLimitSec) = true
MOI.supports(::Optimizer, ::MOI.ObjectiveSense) = true
MOI.supports(::Optimizer, ::MOI.ObjectiveFunction{F}) where {F <: Union{VI, SAF, SQF}} = true

function MOI.set(m::Optimizer, ::MOI.NLPBlock, nlp_data::MOI.NLPBlockData)
    if nlp_data.has_objective
        m._input_problem._objective = nothing
    end
    m._input_problem._nlp_data = nlp_data
end

function MOI.set(m::Optimizer, ::MOI.ObjectiveFunction{T}, f::T) where T <: Union{VI,SAF,SQF}
    check_inbounds!(m, f)
    m._input_problem._objective = f
end
MOI.get(m::Optimizer, ::MOI.ObjectiveFunction{T}) where T <: Union{VI,SAF,SQF} = m._input_problem._objective
MOI.get(m::Optimizer, ::MOI.ObjectiveFunctionType) = typeof(m._input_problem._objective)

MOI.set(m::Optimizer, ::MOI.ObjectiveSense, s::MOI.OptimizationSense) = m._input_problem._optimization_sense = s
MOI.get(m::Optimizer, ::MOI.ObjectiveSense) = m._input_problem._optimization_sense