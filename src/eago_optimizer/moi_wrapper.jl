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

##### Utilities for checking that JuMP model contains variables used in expression
function check_inbounds!(m::Optimizer, vi::VI)
    if !(1 <= vi.value <= m._input_problem._variable_count)
        error("Invalid variable index $vi. ($(m._input_problem._variable_count) variables in the model.)")
    end
    return
end

function check_inbounds!(m::Optimizer, aff::SAF)
    for term in aff.terms
        check_inbounds!(m, term.variable)
    end
    return
end

function check_inbounds!(m::Optimizer, quad::SQF)
    for term in quad.affine_terms
        check_inbounds!(m, term.variable)
    end
    for term in quad.quadratic_terms
        check_inbounds!(m, term.variable_1)
        check_inbounds!(m, term.variable_2)
    end
    return
end

function check_inbounds!(m::Optimizer, vov::VECOFVAR)
    for vi in vov.variables
        check_inbounds!(m, vi)
    end
    return
end

##### Access variable information from MOI variable index
has_upper_bound(m::Optimizer, vi::MOI.VariableIndex) = m._input_problem._variable_info[vi.value].has_upper_bound
has_lower_bound(m::Optimizer, vi::MOI.VariableIndex) = m._input_problem._variable_info[vi.value].has_lower_bound
is_fixed(m::Optimizer, vi::MOI.VariableIndex) = m._input_problem._variable_info[vi.value].is_fixed
is_integer(m::Optimizer, i::Int) = is_integer(m._input_problem._variable_info[i])

##### Add unconstrained variables
function MOI.add_variable(m::Optimizer)
    m._input_problem._variable_count += 1
    push!(m._input_problem._variable_info, VariableInfo{Float64}())
    return VI(m._input_problem._variable_count)
end

##### Supports function and add_constraint for single variable functions
const VAR_SETS = Union{LT, GT, ET, ZO, MOI.Integer}
MOI.supports_constraint(::Optimizer, ::Type{VI}, ::Type{S}) where {S <: VAR_SETS} = true

function MOI.add_constraint(m::Optimizer, v::VI, s::T) where T <: VAR_SETS
    v = v.variable
    check_inbounds!(m, v)
    vi = m._input_problem._variable_info[v.value]
    m._input_problem._variable_info[v.value] = VariableInfo(vi, s)
    return CI{VI, T}(v.value)
end

##### Supports function and add_constraint for scalar affine functions
const INEQ_SETS = Union{LT, GT, ET}
MOI.supports_constraint(::Optimizer, ::Type{SAF}, ::Type{S}) where {S <: INEQ_SETS} = true

macro define_addconstraint_linear(function_type, set_type, array_name)
    quote
        function MOI.add_constraint(m::Optimizer, func::$function_type, set::$set_type)
            check_inbounds!(m, func)
            m._input_problem._last_constraint_index += 1
            civ = m._input_problem._last_constraint_index
            push!(m._input_problem.$(array_name), (func, set, civ))
            indx = CI{$function_type, $set_type}(civ)
            return indx
        end
    end
end
@define_addconstraint_linear SAF LT _linear_leq_constraints
@define_addconstraint_linear SAF GT _linear_geq_constraints
@define_addconstraint_linear SAF ET _linear_eq_constraints

##### Supports function and add_constraint for scalar quadratic functions
MOI.supports_constraint(::Optimizer, ::Type{SQF}, ::Type{S}) where {S <: INEQ_SETS} = true

macro define_addconstraint_quadratic(function_type, set_type, array_name)
    quote
        function MOI.add_constraint(m::Optimizer, func::$function_type, set::$set_type)
            check_inbounds!(m, func)
            m._input_problem._last_constraint_index += 1
            civ = m._input_problem._last_constraint_index
            push!(m._input_problem.$(array_name), (func, set, civ))
            indx = CI{$function_type, $set_type}(civ)
            return indx
        end
    end
end
@define_addconstraint_quadratic SQF LT _quadratic_leq_constraints
@define_addconstraint_quadratic SQF GT _quadratic_geq_constraints
@define_addconstraint_quadratic SQF ET _quadratic_eq_constraints

##### Supports function and add_constraint for conic functions
#=
const CONE_SETS = Union{SOC}
MOI.supports_constraint(::Optimizer, ::Type{VECOFVAR}, ::Type{S}) where {S <: CONE_SETS} = true

function MOI.add_constraint(m::Optimizer, func::VECOFVAR, set::SOC)

    if length(func.variables) !== set.dimension
        error("Dimension of $(s) does not match number of terms in $(f)")
    end

    check_inbounds!(m, func)
    push!(m._input_problem._conic_second_order, (func, set))
    m._input_problem._last_constraint_index += 1
    m._input_problem._conic_second_order_count += 1

    return CI{VECOFVAR, SOC}(m._input_problem._last_constraint_index)
end
=#

function MOI.get(m::Optimizer{R,S,T}, v::MOI.ConstraintPrimal, c::CI{VI,<:Any}) where {R,S,T}
    return MOI.get(m, MOI.VariablePrimal(), MOI.VariableIndex(c.value))
end

function MOI.get(m::Optimizer{R,S,T}, v::MOI.ConstraintPrimal, c::Union{CI{SAF,Q},CI{SQF,Q}}) where {R,S,T,Q}
    return m._global_optimizer._constraint_primal[c.value]
end

function MOI.empty!(m::Optimizer{R,S,T}) where {R,S,T}

    MOI.empty!(m.subsolver_block)
    MOI.empty!(m._global_optimizer)
    m._input_problem = InputProblem()
    m._working_problem = ParsedProblem()

    m._termination_status_code = MOI.OPTIMIZE_NOT_CALLED
    m._result_status_code = MOI.OTHER_RESULT_STATUS

    # set constructor reset on empty! and  to zero in initial parse! in parse.jl
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

    flag &= m._termination_status_code === MOI.OPTIMIZE_NOT_CALLED
    flag &= m._result_status_code === MOI.OTHER_RESULT_STATUS

    # set constructor reset on empty! and  to zero in initial parse! in parse.jl
    flag &= m._run_time == 0.0

    flag &= m._objective_value  == -Inf
    flag &= m._objective_bound  ==  Inf
    flag &= m. _relative_gap    == Inf
    flag &= m._iteration_count  == 0
    flag &= m._node_count       == 0

    return flag
end

MOI.copy_to(model::Optimizer, src::MOI.ModelLike) = MOIU.default_copy_to(model, src)

#####
#####
##### Set & get attributes of model
#####
#####
MOI.supports(::Optimizer, ::MOI.Silent) = true

function MOI.set(m::Optimizer, ::MOI.Silent, value)
    if value
        m._parameters.verbosity = 0
        m._parameters.log_on = false
    else
        m._parameters.verbosity = 1
    end
    return
end

function MOI.set(m::Optimizer, ::MOI.TimeLimitSec, value::Nothing)
    m._parameters.time_limit = Inf
    return
end

function MOI.set(m::Optimizer, ::MOI.TimeLimitSec, value::Float64)
    m._parameters.time_limit = value
    return
end

MOI.get(m::Optimizer, ::MOI.Silent) = m._parameters.verbosity == 0

function MOI.get(m::Optimizer, ::MOI.ListOfVariableIndices)
    return [MOI.VariableIndex(i) for i = 1:length(m._input_problem._variable_info)]
end

MOI.get(m::Optimizer, ::MOI.ObjectiveValue) = m._objective_value

MOI.get(m::Optimizer, ::MOI.NumberOfVariables) = m._input_problem._variable_count
MOI.get(m::Optimizer, ::MOI.ObjectiveBound) = m._objective_bound
function MOI.get(m::Optimizer, ::MOI.RelativeGap)
    b = MOI.get(m, MOI.ObjectiveBound())
    v = MOI.get(m, MOI.ObjectiveValue())
    return relative_gap(b,v)
end

MOI.get(m::Optimizer, ::MOI.SolverName) = "EAGO: Easy Advanced Global Optimization"
MOI.get(m::Optimizer, ::MOI.TerminationStatus) = m._termination_status_code
MOI.get(m::Optimizer, ::MOI.PrimalStatus) = m._result_status_code
MOI.get(m::Optimizer, ::MOI.SolveTimeSec) = m._run_time
MOI.get(m::Optimizer, ::MOI.NodeCount) = m._node_count
MOI.get(m::Optimizer, ::MOI.ResultCount) = (m._result_status_code === MOI.FEASIBLE_POINT) ? 1 : 0
MOI.get(m::Optimizer, ::MOI.TimeLimitSec) = m._parameters.time_limit

function MOI.get(model::Optimizer, ::MOI.VariablePrimal, vi::MOI.VariableIndex)
    check_inbounds!(model, vi)
    return model._global_optimizer._continuous_solution[vi.value]
end

const EAGO_PARAMETERS = fieldnames(EAGOParameters)

_to_sym(d) = error("EAGO only supports raw parameters with Symbol or String names.")
_to_sym(d::String) = Symbol(d)
_to_sym(d::Symbol) = d

function MOI.get(m::Optimizer, p::MOI.RawOptimizerAttribute )
    s = _to_sym(p.name)
    s in EAGO_PARAMETERS ? getfield(m._parameters, s) : getfield(m, s)
end
function MOI.set(m::Optimizer, p::MOI.RawOptimizerAttribute , x)
    s = _to_sym(p.name)
    if (s == :relaxed_optimizer) || (s == :upper_optimizer)
        setfield!(m, s, Incremental(x))
    else
        s in EAGO_PARAMETERS ? setfield!(m._parameters, s, x) : setfield!(m, s, x)
    end
    return
end


#####
#####
##### Support, set, and evaluate objective functions
#####
#####
MOI.supports(::Optimizer, ::MOI.TimeLimitSec) = true
MOI.supports(::Optimizer, ::MOI.ObjectiveSense) = true
MOI.supports(::Optimizer, ::MOI.ObjectiveFunction{F}) where {F <: Union{VI, SAF, SQF}} = true

function MOI.set(m::Optimizer, ::MOI.NLPBlock, nlp_data::MOI.NLPBlockData)
    if nlp_data.has_objective
        m._input_problem._objective = nothing
    end
    m._input_problem._nlp_data = nlp_data
    return
end

function MOI.set(m::Optimizer, ::MOI.ObjectiveFunction{T}, f::T) where T <: Union{VI,SAF,SQF}
    check_inbounds!(m, f)
    m._input_problem._objective = f
    return
end

function MOI.set(m::Optimizer, ::MOI.ObjectiveSense, s::MOI.OptimizationSense)
    m._input_problem._optimization_sense = s
    return
end
