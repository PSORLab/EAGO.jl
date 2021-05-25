# Copyright (c) 2018: Matthew Wilhelm & Matthew Stuber.
# This code is licensed under MIT license (see LICENSE.md for full details)
#############################################################################
# EAGO
# A development environment for robust and global optimization
# See https://github.com/PSORLab/EAGO.jl
#############################################################################
# src/eago_optimizer/constraints.jl
# Defines constraints supported Optimizer and how to store them.
#############################################################################

##### Access variable information from MOI variable index
has_upper_bound(m::Optimizer, vi::MOI.VariableIndex) = m._input_problem._variable_info[vi.value].has_upper_bound
has_lower_bound(m::Optimizer, vi::MOI.VariableIndex) = m._input_problem._variable_info[vi.value].has_lower_bound
is_fixed(m::Optimizer, vi::MOI.VariableIndex) = m._input_problem._variable_info[vi.value].is_fixed
is_integer(m::Optimizer, i::Int64) = is_integer(m._input_problem._variable_info[i])

##### Add unconstrained variables
function MOI.add_variable(m::Optimizer)
    m._input_problem._variable_count += 1
    push!(m._input_problem._variable_info, VariableInfo())
    return VI(m._input_problem._variable_count)
end
MOI.add_variables(m::Optimizer, n::Int) = [MOI.add_variable(m) for i in 1:n]

##### Supports function and add_constraint for single variable functions
const INEQ_SETS = Union{LT, GT, ET}
MOI.supports_constraint(::Optimizer, ::Type{SV}, ::Type{S}) where {S <: INEQ_SETS} = true

function MOI.add_constraint(m::Optimizer, v::SV, s::T) where T <: INEQ_SETS
    v = v.variable
    check_inbounds!(m, v)
    vi = m._input_problem._variable_info[v.value]
    m._input_problem._variable_info[v.value] = VariableInfo(vi, s)
    return CI{SV, T}(v.value)
end

##### Supports function and add_constraint for scalar affine functions
MOI.supports_constraint(::Optimizer, ::Type{SAF}, ::Type{S}) where {S <: INEQ_SETS} = true

macro define_addconstraint_linear(function_type, set_type, array_name, count_name)
    quote
        function MOI.add_constraint(m::Optimizer, func::$function_type, set::$set_type)
            check_inbounds!(m, func)
            push!(m._input_problem.$(array_name), (func, set))
            m._input_problem._last_constraint_index += 1
            m._input_problem.$(count_name) += 1
            indx = CI{$function_type, $set_type}(m._input_problem._last_constraint_index)
            return indx
        end
    end
end

@define_addconstraint_linear SAF LT _linear_leq_constraints _linear_leq_count
@define_addconstraint_linear SAF GT _linear_geq_constraints _linear_geq_count
@define_addconstraint_linear SAF ET _linear_eq_constraints _linear_eq_count

##### Supports function and add_constraint for scalar quadratic functions
MOI.supports_constraint(::Optimizer, ::Type{SQF}, ::Type{S}) where {S <: INEQ_SETS} = true

macro define_addconstraint_quadratic(function_type, set_type, array_name, count_name)
    quote
        function MOI.add_constraint(m::Optimizer, func::$function_type, set::$set_type)
            check_inbounds!(m, func)
            push!(m._input_problem.$(array_name), (func, set))
            m._input_problem._last_constraint_index += 1
            m._input_problem.$(count_name) += 1
            indx = CI{$function_type, $set_type}(m._input_problem._last_constraint_index)
            return indx
        end
    end
end

@define_addconstraint_quadratic SQF LT _quadratic_leq_constraints _quadratic_leq_count
@define_addconstraint_quadratic SQF GT _quadratic_geq_constraints _quadratic_geq_count
@define_addconstraint_quadratic SQF ET _quadratic_eq_constraints _quadratic_eq_count

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
