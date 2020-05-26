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
# src/eago_optimizer/constraints.jl
# Defines constraints supported Optimizer and how to store them.
#############################################################################


##### Supports function and add_constraint for scalar affine functions
MOI.supports_constraint(::Optimizer, ::Type{SAF}, ::Type{S}) where {S <: INEQ_SETS} = true

macro define_addconstraint_linear(function_type, set_type, array_name, count_name)
    quote
        function MOI.add_constraint(m::Optimizer, func::$function_type, set::$set_type)
            check_inbounds!(m, func)
            push!(m.$(array_name), (func, set))
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
            for i in func.affine_terms m.branch_variable[i.variable_index.value] = true end
            for i in func.quadratic_terms
                m.branch_variable[i.variable_index_1.value] = true
                m.branch_variable[i.variable_index_2.value] = true
            end
            push!(m.$(array_name), (func, set))
            m._input_problem._last_constraint_index += 1
            m._input_problem.$(count_name) += 1
            indx = CI{$function_type, $set_type}(m._input_problem._last_constraint_index)
            return indx
        end
    end
end

@define_addconstraint_quadratic SQF LT _quadratic_leq_constraints _quadratic_leq_constraints
@define_addconstraint_quadratic SQF GT _quadratic_geq_constraints _quadratic_geq_constraints
@define_addconstraint_quadratic SQF ET _quadratic_eq_constraints _quadratic_eq_constraints

##### Supports function and add_constraint for conic functions
const CONE_SETS = Union{MOI.NormInfinityCone, MOI.NormOneCone, MOI.SecondOrderCone, MOI.RotatedSecondOrderCone,
                        MOI.GeometricMeanCone, MOI.ExponentialCone, MOI.DualExponentialCone, MOI.PowerCone,
                        MOI.DualPowerCone, MOI.RelativeEntropyCone, MOI.NormSpectralCone, MOI.NormNuclearCone}
MOI.supports_constraint(::Optimizer, ::Type{VECOFVAR}, ::Type{S}) where {S <: CONE_SETS} = true

macro define_addconstraint_cone(set_type, array_name, count_name)
    quote
        function MOI.add_constraint(m::Optimizer, func::VECOFVAR, set::$set_type)
            if length(func.variables) !== dimension(set)
                error("Dimension of $(s) does not match number of terms in $(f)")
            end
            check_inbounds!(m, func)
            push!(m.$(array_name), (func, set))
            m._input_problem._last_constraint_index += 1
            m._input_problem.$(count_name) += 1
            return CI{VECOFVAR, $set_type}(m._input_problem._last_constraint_index)
        end
    end
end

@define_addconstraint_cone MOI.SecondOrderCone _conic_second_order _conic_second_order_count
