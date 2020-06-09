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
# src/eeago_optimizer/optimize/optimize_convex.jl
# Contains the solve_local_nlp! routine which computes the optimal value
# of a convex function. This is used to compute the upper bound in the
# branch and bound routine. A number of utility function required for
# solve_local_nlp! are also included.
#############################################################################

"""
$(SIGNATURES)

Shifts the resulting local nlp objective value `f*` by `(1.0 + relative_tolerance/100.0)*f* + absolute_tolerance/100.0`.
This assumes that the local solvers relative tolerance and absolute tolerance is significantly lower than the global
tolerance (local problem is minimum).
"""
function stored_adjusted_upper_bound!(d::Optimizer, v::Float64)
    adj_atol = d._parameters.absolute_tolerance/100.0
    adj_rtol = d._parameters.relative_tolerance/100.0
    if v > 0.0
        d._upper_objective_value = v*(1.0 + adj_rtol) + adj_atol
    else
        d._upper_objective_value = v*(1.0 - adj_rtol) + adj_atol
    end

    return nothing
end

# translates quadratic cone
function add_soc_constraints_as_quad!(m::Optimizer, opt::T) where T

    for (func, set) in m._input_problem._conic_second_order
        # quadratic cone implies variable[1] >= 0.0, bounds contracted accordingly in initial_parse!
        quad_terms = SQT[SQT((), func.variables[i], func.variables[i]) for i = 1:length(func.variables)]
        sqf = SQF(SQT[], SAF[], 0.0)
        MOI.add_constraint(opt, sqf, LT_ZERO)
    end

    return nothing
end

"""
$(SIGNATURES)

Constructs and solves the problem locally on on node `y` updated the upper
solution informaton in the optimizer.
"""
function solve_local_nlp!(m::Optimizer)

    upper_optimizer = m.upper_optimizer
    MOI.empty!(upper_optimizer)

    upper_variables = m._upper_variables
    for i = 1:m._working_problem._variable_count
        @inbounds upper_variables[i] = MOI.add_variable(upper_optimizer)
    end

    n = m._current_node
    sol_to_branch_map = m._sol_to_branch_map
    lower_variable_bounds = n.lower_variable_bounds
    upper_variable_bounds = n.upper_variable_bounds
    variable_info = m._input_problem._variable_info

    lvb = 0.0
    uvb = 0.0
    x0 = 0.0

    for i = 1:m._input_problem._variable_count
        vinfo = @inbounds variable_info[i]
        single_variable = MOI.SingleVariable(@inbounds upper_variables[i])

        if vinfo.branch_on === BRANCH
            if vinfo.is_integer
            else
                indx = @inbounds sol_to_branch_map[i]
                lvb  = @inbounds lower_variable_bounds[indx]
                uvb  = @inbounds upper_variable_bounds[indx]
                if vinfo.is_fixed
                    MOI.add_constraint(upper_optimizer, single_variable, ET(lvb))

                elseif vinfo.has_lower_bound
                    if vinfo.has_upper_bound
                        MOI.add_constraint(upper_optimizer, single_variable, LT(uvb))
                        MOI.add_constraint(upper_optimizer, single_variable, GT(lvb))

                    else
                        MOI.add_constraint(upper_optimizer, single_variable, GT(lvb))

                    end
                elseif vinfo.has_upper_bound
                    MOI.add_constraint(upper_optimizer, single_variable, LT(uvb))

                end
            end
            x0 = 0.5*(lvb + uvb)
            upper_variable_index = @inbounds upper_variables[i]
            MOI.set(upper_optimizer, MOI.VariablePrimalStart(), upper_variable_index, x0)

        else
            # not branch variable
            if vinfo.is_integer
            else
                lvb  = vinfo.lower_bound
                uvb  = vinfo.upper_bound
                if vinfo.is_fixed
                    MOI.add_constraint(upper_optimizer, single_variable, ET(lvb))

                elseif vinfo.has_lower_bound
                    if vinfo.has_upper_bound
                        MOI.add_constraint(upper_optimizer, single_variable, LT(uvb))
                        MOI.add_constraint(upper_optimizer, single_variable, GT(lvb))

                    else
                        MOI.add_constraint(upper_optimizer, single_variable, GT(lvb))

                    end
                elseif vinfo.has_upper_bound
                    MOI.add_constraint(upper_optimizer, single_variable, LT(uvb))
                end
                x0 = 0.5*(lvb + uvb)
                upper_variable_index = @inbounds upper_variables[i]
                MOI.set(upper_optimizer, MOI.VariablePrimalStart(), upper_variable_index, x0)
            end
        end
    end

    # Add linear and quadratic constraints to model
    add_linear_constraints!(m, upper_optimizer)

    for (func, set) in m._input_problem._quadratic_leq_constraints
        MOI.add_constraint(upper_optimizer, func, set)
    end
    for (func, set) in m._input_problem._quadratic_geq_constraints
        MOI.add_constraint(upper_optimizer, func, set)
    end
    for (func, set) in m._input_problem._quadratic_eq_constraints
        MOI.add_constraint(upper_optimizer, func, set)
    end

    if MOI.supports_constraint(upper_optimizer, VECOFVAR, SOC)
        add_soc_constraints!(m, upper_optimizer)
    else
        add_soc_constraints_as_quad!(m, upper_optimizer)
    end

    # Add nonlinear evaluation block
    MOI.set(upper_optimizer, MOI.NLPBlock(), m._working_problem._nlp_data)
    MOI.set(upper_optimizer, MOI.ObjectiveSense(), MOI.MIN_SENSE)

    # set objective as NECESSARY
    add_sv_or_aff_obj!(m, upper_optimizer)
    if m._input_problem._objective_type === SCALAR_QUADRATIC
        MOI.set(upper_optimizer, MOI.ObjectiveFunction{SQF}(), m._input_problem._objective_sqf)
    end

    # Optimizes the object
    MOI.optimize!(upper_optimizer)

    # Process output info and save to CurrentUpperInfo object
    m._upper_termination_status = MOI.get(upper_optimizer, MOI.TerminationStatus())
    m._upper_result_status = MOI.get(upper_optimizer, MOI.PrimalStatus())

    if is_feasible_solution(m._upper_termination_status, m._upper_result_status)
        m._upper_feasibility = true
        value = MOI.get(upper_optimizer, MOI.ObjectiveValue())
        stored_adjusted_upper_bound!(m, value)
        m._best_upper_value = min(value, m._best_upper_value)
        m._upper_solution .= MOI.get(upper_optimizer, MOI.VariablePrimal(), upper_variables)

    else
        m._upper_feasibility = false
        m._upper_objective_value = Inf

    end

    return nothing
end

function optimize!(::Val{DIFF_CVX}, m::Optimizer)

    solve_local_nlp!(m)

    return nothing
end
