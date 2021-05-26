"""
$(SIGNATURES)

Constructs and solves the relaxation using the default EAGO relaxation scheme
and optimizer on node `y`.
"""
function lower_problem!(t::ExtensionType, m::Optimizer)

    n = m._current_node
    wp = m._working_problem

    wp._relaxed_evaluator.is_post = m._parameters.subgrad_tighten
    if !m._obbt_performed_flag
        if m._nonlinear_evaluator_created
            set_node!(wp._relaxed_evaluator, n)
            set_reference_point!(m)
            fill!(wp._relaxed_evaluator.subexpressions_eval, false)
        end
        update_relaxed_problem_box!(m)
    end
    _set_has_value!(wp._objective_nl, false)
    wp._relaxed_evaluator.interval_intersect = false

    if !m._obbt_performed_flag
        relax_constraints!(m, 1)
    end
    relax_objective!(m, 1)

    # Optimizes the object
    relaxed_optimizer = m.relaxed_optimizer

    MOI.set(relaxed_optimizer, MOI.ObjectiveSense(), MOI.MIN_SENSE)

    MOI.optimize!(relaxed_optimizer)

    m._lower_termination_status = MOI.get(relaxed_optimizer, MOI.TerminationStatus())
    m._lower_primal_status = MOI.get(relaxed_optimizer, MOI.PrimalStatus())
    m._lower_dual_status = MOI.get(relaxed_optimizer, MOI.DualStatus())
    status = relaxed_problem_status(m._lower_termination_status,
                                    m._lower_primal_status,
                                    m._lower_dual_status)

    if status == RRS_OPTIMAL
        set_dual!(m)
        m._cut_add_flag = true
        m._lower_feasibility = true
        m._lower_objective_value = MOI.get(relaxed_optimizer, MOI.ObjectiveValue())
        for i = 1:m._working_problem._variable_count
             m._lower_solution[i] = MOI.get(relaxed_optimizer, MOI.VariablePrimal(), m._relaxed_variable_index[i])
        end
    elseif status == RRS_INFEASIBLE
        m._cut_add_flag = false
        m._lower_feasibility  = false
        m._lower_objective_value = -Inf
    #elseif status == RRS_DUAL_FEASIBLE
    else
        fallback_interval_lower_bound!(m, n)
    end

    return nothing
end
lower_problem!(m::Optimizer) = lower_problem!(m.ext_type, m)

"""
$(SIGNATURES)

Updates the internal storage in the optimizer after a valid feasible cut is added.
"""
function cut_update!(m::Optimizer)

    m._cut_feasibility = true

    relaxed_optimizer = m.relaxed_optimizer
    obj_val = MOI.get(relaxed_optimizer, MOI.ObjectiveValue())
    prior_obj_val = (m._cut_iterations == 2) ? m._lower_objective_value : m._cut_objective_value

    m._cut_add_flag = true
    m._lower_termination_status = m._cut_termination_status
    m._lower_primal_status = m._cut_primal_status
    m._lower_dual_status = m._cut_dual_status
    m._cut_solution[:] = MOI.get(relaxed_optimizer, MOI.VariablePrimal(), m._relaxed_variable_index)

    if prior_obj_val < obj_val
        m._cut_objective_value = obj_val
        m._lower_objective_value = obj_val
        set_dual!(m)
        copyto!(m._lower_solution, m._cut_solution)

    else
        m._cut_objective_value = prior_obj_val
        m._lower_objective_value = prior_obj_val
        m._cut_add_flag = false
    end

    return nothing
end

"""
$(SIGNATURES)

Checks if a cut should be added and computes a new reference point to add the
cut at. If no cut should be added the constraints not modified in place are
deleted from the relaxed optimizer and the solution is compared with the
interval lower bound. The best lower bound is then used.
"""
function cut_condition(t::ExtensionType, m::Optimizer)

    # always add cut if below the minimum iteration limit, otherwise add cut
    # the number of cuts is less than the maximum and the distance between
    # prior solutions exceeded a tolerance.
    continue_cut_flag = m._cut_add_flag
    continue_cut_flag &= (m._cut_iterations < m._parameters.cut_max_iterations)

    # compute distance between prior solutions and compare to tolerances
    n = m._current_node
    ns_indx = m._branch_to_sol_map

    cvx_factor =  m._parameters.cut_cvx
    xsol = (m._cut_iterations > 1) ? m._cut_solution[ns_indx] : m._lower_solution[ns_indx]
    xnew = (1.0 - cvx_factor)*mid(n) + cvx_factor*xsol

    continue_cut_flag &= (norm((xsol - xnew)/diam(n), 1) > m._parameters.cut_tolerance)
    continue_cut_flag |= (m._cut_iterations < m._parameters.cut_min_iterations)

    # update reference point for new cut
    if continue_cut_flag
        copyto!(m._current_xref, xnew)
        if m._nonlinear_evaluator_created
            set_reference_point!(m)
            fill!(m._working_problem._relaxed_evaluator.subexpressions_eval, false)
        end
    end

    # check to see if interval bound is preferable and replaces the objective
    # value with the interval value if so. Any available dual values are then
    # set to zero since the interval bounds are by definition constant
    if m._lower_feasibility && !continue_cut_flag
        objective_lo = -Inf
        obj_type = m._working_problem._objective_type
        if obj_type === SINGLE_VARIABLE
            var_index = m._working_problem._objective_sv.variable.value
            if m._branch_variables[var_index]
                obj_indx = m._sol_to_branch_map[var_index]
                lower_variable_bnd = n.lower_variable_bounds[obj_indx]
                if !isinf(lower_variable_bnd)
                    objective_lo = lower_variable_bnd
                end
            end

        elseif obj_type === SCALAR_AFFINE
            objective_lo = lower_interval_bound(m, m._working_problem._objective_saf_parsed, n)

        elseif obj_type === SCALAR_QUADRATIC
            objective_lo = lower_interval_bound(m, m._working_problem._objective_sqf, n)

        elseif obj_type === NONLINEAR
            objective_lo = lower_interval_bound(m, m._working_problem._objective_nl, n)

        end

        if objective_lo > m._lower_objective_value
            m._lower_objective_value = objective_lo
            fill!(m._lower_lvd, 0.0)
            fill!(m._lower_uvd, 0.0)
        end
    end

    m._cut_iterations += 1

    return continue_cut_flag
end
cut_condition(m::Optimizer) = cut_condition(m.ext_type, m)

"""
$(SIGNATURES)

Adds a cut for each constraint and the objective function to the subproblem.
"""
function add_cut!(t::ExtensionType, m::Optimizer)

    fill!(m._working_problem._relaxed_evaluator.subexpressions_eval, false)
    m._working_problem._relaxed_evaluator.is_first_eval = true
    m._working_problem._relaxed_evaluator.is_intersect = false
    m._new_eval_objective = true
    m._new_eval_constraint = true

    relax_constraints!(m, m._cut_iterations)
    relax_objective!(m, m._cut_iterations)

    # Optimizes the object
    relaxed_optimizer = m.relaxed_optimizer
    MOI.optimize!(relaxed_optimizer)

    m._cut_termination_status = MOI.get(relaxed_optimizer, MOI.TerminationStatus())
    m._cut_primal_status = MOI.get(relaxed_optimizer, MOI.PrimalStatus())
    m._cut_dual_status = MOI.get(relaxed_optimizer, MOI.DualStatus())
    status = relaxed_problem_status(m._cut_termination_status,
                                    m._cut_primal_status,
                                    m._cut_dual_status)

    if status == RRS_OPTIMAL
        cut_update!(m)

    elseif valid_flag
        m._cut_add_flag = false
        m._lower_feasibility  = false
        m._lower_objective_value = -Inf

    else
        m._cut_add_flag = false
    end

    return nothing
end
add_cut!(m::Optimizer) = add_cut!(m.ext_type, m)
