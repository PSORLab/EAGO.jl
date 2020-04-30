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

macro define_addconstraint_linear(function_type, set_type, array_name)
    quote
        function MOI.add_constraint(m::Optimizer, func::$function_type, set::$set_type)
            check_inbounds!(m, func)
            push!(m.$(array_name), (func, set))
            m._last_constraint_index += 1
            indx = CI{$function_type, $set_type}(m._last_constraint_index)
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
            for i in func.affine_terms m.branch_variable[i.variable_index.value] = true end
            for i in func.quadratic_terms
                m.branch_variable[i.variable_index_1.value] = true
                m.branch_variable[i.variable_index_2.value] = true
            end
            push!(m.$(array_name), (func, set))
            m._last_constraint_index += 1
            indx = CI{$function_type, $set_type}(m._last_constraint_index)
            return indx
        end
    end
end

@define_addconstraint_quadratic SQF LT _quadratic_leq_constraints
@define_addconstraint_quadratic SQF GT _quadratic_geq_constraints
@define_addconstraint_quadratic SQF ET _quadratic_eq_constraints


##### Supports function and add_constraint for conic functions
const CONE_SETS = Union{MOI.NormInfinityCone, MOI.NormOneCone, MOI.SecondOrderCone, MOI.RotatedSecondOrderCone,
                        MOI.GeometricMeanCone, MOI.ExponentialCone, MOI.DualExponentialCone, MOI.PowerCone,
                        MOI.DualPowerCone, MOI.RelativeEntropyCone, MOI.NormSpectralCone, MOI.NormNuclearCone}
MOI.supports_constraint(::Optimizer, ::Type{VECOFVAR}, ::Type{S}) where {S <: CONE_SETS} = true

macro define_addconstraint_cone(set_type, array_name)
    quote
        function MOI.add_constraint(m::Optimizer, func::VECOFVAR, set::$set_type)
            if length(func.variables) != dimension(set)
                error("Dimension of $(s) does not match number of terms in $(f)")
            end
            check_inbounds!(m, func)
            push!(m.$(array_name), (func, set))
            m._last_constraint_index += 1
            return CI{VECOFVAR, $set_type}(m._last_constraint_index)
        end
    end
end

@define_addconstraint_cone MOI.NormInfinityCone _conic_norm_infinity
@define_addconstraint_cone MOI.NormOneCone _conic_norm_one
@define_addconstraint_cone MOI.SecondOrderCone _conic_second_order
@define_addconstraint_cone MOI.RotatedSecondOrderCone _conic_rotated_second_order
@define_addconstraint_cone MOI.GeometricMeanCone _conic_geometric_mean
@define_addconstraint_cone MOI.ExponentialCone _conic_exponential
@define_addconstraint_cone MOI.DualExponentialCone _conic_dual_exponential
@define_addconstraint_cone MOI.PowerCone _conic_power_cone
@define_addconstraint_cone MOI.DualPowerCone _conic_dual_power
@define_addconstraint_cone MOI.RelativeEntropyCone _conic_relative_entropy
@define_addconstraint_cone MOI.NormSpectralCone _conic_norm_spectral
@define_addconstraint_cone MOI.NormNuclearCone _conic_norm_nuclear
