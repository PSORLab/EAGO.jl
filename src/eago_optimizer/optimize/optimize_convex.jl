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
    nothing
end

"""
$(SIGNATURES)

Constructs and solves the problem locally on on node `y` updated the upper
solution informaton in the optimizer.
"""
function single_nlp_solve!(m::Optimizer)

    upper_optimizer = m.upper_optimizer
    MOI.empty!(upper_optimizer)

    upper_variables = MOI.add_variables(upper_optimizer, m._variable_number)

    n = m._current_node
    lower_variable_bounds = n.lower_variable_bounds
    upper_variable_bounds = n.upper_variable_bounds
    variable_info = m._input_problem._variable_info

    lvb = 0.0
    uvb = 0.0
    x0 = 0.0
    for i = 1:m._input_problem._variable_number
        vinfo = @inbounds variable_info[i]
        single_variable = MOI.SingleVariable(@inbounds upper_variables[i])
        if vinfo.is_integer
        else
            lvb = @inbounds lower_variable_bounds[i]
            uvb = @inbounds upper_variable_bounds[i]
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
            MOI.set(upper_optimizer, MOI.VariablePrimalStart(), @inbounds upper_variables[i], x0)
        end
    end

    # Add linear and quadratic constraints to model
    for (func, set) in m._input_problem._linear_leq_constraints
         MOI.add_constraint(upper_optimizer, func, set)
    end
    for (func, set) in m._input_problem._linear_geq_constraints
        MOI.add_constraint(upper_optimizer, func, set)
    end
    for (func, set) in m._input_problem._linear_eq_constraints
        MOI.add_constraint(upper_optimizer, func, set)
    end

    for (func, set) in m._input_problem._quadratic_leq_constraints
        MOI.add_constraint(upper_optimizer, func, set)
    end
    for (func, set) in m._input_problem._quadratic_geq_constraints
        MOI.add_constraint(upper_optimizer, func, set)
    end
    for (func, set) in m._input_problem._quadratic_eq_constraints
        MOI.add_constraint(upper_optimizer, func, set)
    end

    # Add nonlinear evaluation block
    MOI.set(upper_optimizer, MOI.NLPBlock(), m._input_problem._nlp_data)
    MOI.set(upper_optimizer, MOI.ObjectiveSense(), MOI.MIN_SENSE)

    if x._input_problem._objective_type === SINGLE_VARIABLE
        MOI.set(upper_optimizer, MOI.ObjectiveFunction{SV}(), m._input_problem._objective_sv)
    elseif x._input_problem._objective_type === SCALAR_AFFINE
        MOI.set(upper_optimizer, MOI.ObjectiveFunction{SAF}(), m._input_problem._objective_saf)
    elseif x._input_problem._objective_type === SCALAR_QUADRATIC
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
        m._upper_solution .= MOI.get(upper_optimizer, MOI.VariablePrimal(), upper_vars)
    else
        m._upper_feasibility = false
        m._upper_objective_value = Inf
    end

    return nothing
end
