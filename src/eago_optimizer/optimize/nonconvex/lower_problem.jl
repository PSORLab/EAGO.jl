
function relax_objective!(m, f::T, check_safe) where T<:Union{SV,SAF}
    MOI.set(m.relaxed_optimizer, MOI.ObjectiveFunction{T}(), f)
    return
end
function relax_objective!(m, f::BufferedQuadraticIneq, check_safe)
    wp = m._working_problem
    finite_cut_generated = affine_relax_quadratic!(f.func, f.buffer, f.saf,
                                                   m._current_node, m._sol_to_branch_map,
                                                   m._current_xref)
    if finite_cut_generated && (!check_safe || is_safe_cut!(m, f.saf))
        copyto!(wp._objective_saf.terms, f.saf.terms)
        wp._objective_saf.constant = f.saf.constant
        MOI.set(m.relaxed_optimizer, MOI.ObjectiveFunction{SAF}(), wp._objective_saf)
    end
    return
end
function relax_objective!(m, f::BufferedNonlinearFunction, check_safe)
    wp = m._working_problem
    new_flag = m._new_eval_objective
    wp._relaxed_evaluator.is_first_eval = new_flag
    finite_cut_generated = affine_relax_nonlinear!(f, wp._relaxed_evaluator, true, new_flag, false)
    wp._relaxed_evaluator.is_first_eval = false

    if finite_cut_generated && (!check_safe || is_safe_cut!(m, f.saf))
        copyto!(wp._objective_saf.terms, f.saf.terms)
        wp._objective_saf.constant = f.saf.constant
        MOI.set(m.relaxed_optimizer, MOI.ObjectiveFunction{SAF}(), wp._objective_saf)
    end
    return
end

"""
$(TYPEDSIGNATURES)

Triggers an evaluation of the objective function and then updates
the affine relaxation of the objective function.
"""
function relax_objective!(t::ExtensionType, m::Optimizer, q::Int)
    check_safe = (q === 1) ? false : m._parameters.cut_safe_on
    relax_objective!(m, m._working_problem._objective, check_safe)
    m._new_eval_objective = false
    return
end
relax_objective!(m::Optimizer, q::Int64) = relax_objective!(m.ext_type, m, q)

function objective_cut!(m::Optimizer, f::SV, check_safe::Bool, ϵ)
    if !isinf(m._global_upper_bound) && (m._objective_cut_ci_sv.value === -1)
        m._objective_cut_ci_sv = CI{SV,LT}(f.variable.value)
    end
    MOI.set(m.relaxed_optimizer, MOI.ConstraintSet(), m._objective_cut_ci_sv, LT(UBD))
    return
end
function objective_cut!(m::Optimizer, f::AffineFunctionIneq, check_safe::Bool, ϵ)
    formulated_constant = wp._objective_saf.constant
    wp._objective_saf.constant = 0.0
    if check_safe && is_safe_cut!(m, wp._objective_saf)
        ci_saf = MOI.add_constraint(m.relaxed_optimizer, wp._objective_saf, LT(UBD - wp._objective_saf.constant + ϵ))
        push!(m._objective_cut_ci_saf, ci_saf)
    end
    wp._objective_saf.constant = formulated_constant
    return
end
function objective_cut!(m::Optimizer, f::BufferedQuadraticIneq, check_safe::Bool)
    finite_cut_generated = affine_relax_quadratic!(f.func, f.buffer,
                                                   f.saf, m._current_node, m._sol_to_branch_map,
                                                   m._current_xref)
    if finite_cut_generated && (!check_safe || is_safe_cut!(m, f.saf))
        copyto!(wp._objective_saf.terms, f.saf.terms)
        wp._objective_saf.constant = 0.0
        ci_saf = MOI.add_constraint(m.relaxed_optimizer, wp._objective_saf, LT(UBD - buffered_sqf.saf.constant + ϵ))
        push!(m._objective_cut_ci_saf, ci_saf)
    end
    return
end
function objective_cut!(m::Optimizer, f::BufferedNonlinearFunction, check_safe::Bool, ϵ)

    UBD = m._global_upper_bound
    relaxed_optimizer = m.relaxed_optimizer
    wp = m._working_problem
    relaxed_evaluator = wp._relaxed_evaluator

    # if the objective cut is the first evaluation of the objective expression
    # then perform a a forward pass
    new_flag = m._new_eval_objective
    relaxed_evaluator.is_first_eval = new_flag
    finite_cut_generated = affine_relax_nonlinear!(f, relaxed_evaluator, true, new_flag, false)

    constraint_tol = m._parameters.absolute_constraint_feas_tolerance
    if finite_cut_generated
        copyto!(wp._objective_saf.terms, f.saf.terms)
        wp._objective_saf.constant = 0.0
        if !check_safe || is_safe_cut!(m,  f.saf)
            # TODO: When we introduce numerically safe McCormick operators we'll need to replace
            # the UBD - buffered_nl.saf.constant with a correctly rounded version. For now,
            # a small factor is added to the UBD calculation initially which should be sufficient.
            ci_saf = MOI.add_constraint(m.relaxed_optimizer, wp._objective_saf, LT(UBD - f.saf.constant + ϵ))
            push!(m._objective_cut_ci_saf, ci_saf)
        end
    end
    return
end

"""
$(FUNCTIONNAME)

Adds linear objective cut constraint to the `x.relaxed_optimizer`.
"""
function objective_cut!(m::Optimizer, check_safe::Bool)
    ϵ = m._parameters.absolute_constraint_feas_tolerance
    if m._parameters.objective_cut_on && m._global_upper_bound < Inf
        objective_cut!(m, m._working_problem._objective, check_safe, ϵ)
        m._new_eval_objective = false
    end
    return
end

"""
    RelaxResultStatus

Status code used internally to determine how to interpret theresults from the
solution of a relaxed problem.
"""
@enum(RelaxResultStatus, RRS_OPTIMAL, RRS_DUAL_FEASIBLE, RRS_INFEASIBLE, RRS_INVALID)

"""
$(SIGNATURES)

Takes an `MOI.TerminationStatusCode` and a `MOI.ResultStatusCode` and returns
the tuple `(valid_result::Bool, feasible::Bool)`. The value `valid_result` is
`true` if the pair of codes prove that either the subproblem solution was solved
to global optimality or the subproblem solution is infeasible. The value of
`feasible` is true if the problem is feasible and false if the problem is infeasible.
"""
function relaxed_problem_status(t::MOI.TerminationStatusCode,
                                p::MOI.ResultStatusCode,
                                d::MOI.ResultStatusCode)

    if (t == MOI.OPTIMAL) && (p == MOI.FEASIBLE_POINT)
        return RRS_OPTIMAL
    elseif t == MOI.INFEASIBLE
        if (p == MOI.INFEASIBILITY_CERTIFICATE) ||
           (p == MOI.NO_SOLUTION) || (p == MOI.UNKNOWN_RESULT_STATUS)
            return RRS_INFEASIBLE
        end
    elseif (t == MOI.INFEASIBLE_OR_UNBOUNDED && p == MOI.NO_SOLUTION)
        return RRS_INFEASIBLE
    end
    (d == MOI.FEASIBLE_POINT) && return RRS_DUAL_FEASIBLE
    return RRS_INVALID
end

"""
$(SIGNATURES)

Updates the relaxed constraint by setting the constraint set of `v == x*`` ,
`xL_i <= x_i`, and `x_i <= xU_i` for each such constraint added to the relaxed
optimizer.
"""
function update_relaxed_problem_box!(m::Optimizer)

    opt = m.relaxed_optimizer
    wp = m._working_problem

    n = m._current_node
    lower_bound = n.lower_variable_bounds
    upper_bound = n.upper_variable_bounds

    relaxed_variable_eq = m._relaxed_variable_eq
    for i = 1:wp._var_eq_count
        constr_indx, node_indx =  relaxed_variable_eq[i]
        MOI.set(opt, MOI.ConstraintSet(), constr_indx, ET( lower_bound[node_indx]))
    end

    relaxed_variable_lt = m._relaxed_variable_lt
    for i = 1:wp._var_leq_count
        constr_indx, node_indx =  relaxed_variable_lt[i]
        MOI.set(opt, MOI.ConstraintSet(), constr_indx, LT( upper_bound[node_indx]))
    end

    relaxed_variable_gt = m._relaxed_variable_gt
    for i = 1:wp._var_geq_count
        constr_indx, node_indx =  relaxed_variable_gt[i]
        MOI.set(opt, MOI.ConstraintSet(), constr_indx, GT( lower_bound[node_indx]))
    end

    return nothing
end

function reset_relaxation!(m::Optimizer)

    m._working_problem._relaxed_evaluator.is_first_eval = true
    fill!(m._working_problem._relaxed_evaluator.subexpressions_eval, false)

    m._new_eval_objective = true
    m._new_eval_constraint = true

    delete_nl_constraints!(m)
    delete_objective_cuts!(m)

    return nothing
end

"""
$(SIGNATURES)

Retrieves the lower and upper duals for variable bounds from the
`relaxed_optimizer` and sets the appropriate values in the
`_lower_lvd` and `_lower_uvd` storage fields.
"""
function set_dual!(m::Optimizer)

    relaxed_optimizer = m.relaxed_optimizer
    relaxed_variable_lt = m._relaxed_variable_lt
    relaxed_variable_gt = m._relaxed_variable_gt

    for i = 1:m._working_problem._var_leq_count
        ci_lt, i_lt = @inbounds relaxed_variable_lt[i]
        @inbounds m._lower_uvd[i_lt] = MOI.get(relaxed_optimizer, MOI.ConstraintDual(), ci_lt)
    end
    for i = 1:m._working_problem._var_geq_count
        ci_gt, i_gt = @inbounds relaxed_variable_gt[i]
        @inbounds m._lower_lvd[i_gt] = MOI.get(relaxed_optimizer, MOI.ConstraintDual(), ci_gt)
    end

    return nothing
end

function interval_objective_bound(m::Optimizer, n::NodeBB)

    interval_objective_bound = bound_objective(m)

    if interval_objective_bound > m._lower_objective_value
        m._lower_objective_value = interval_objective_bound
        fill!(m._lower_lvd, 0.0)
        fill!(m._lower_uvd, 0.0)
        m._cut_add_flag = false
        return true
    end

    return false
end

"""
$(SIGNATURES)

A fallback lower bounding problem that consists of an natural interval extension
calculation. This is called when the optimizer used to compute the lower bound
does not return a termination and primal status code indicating that it
successfully solved the relaxation to a globally optimal point.
"""
function fallback_interval_lower_bound!(m::Optimizer, n::NodeBB)

    feasible_flag = true
    wp = m._working_problem

    if !cp_condition(m)
        feas = feas && all(f -> is_feasible(m, f, n), wp._saf_leq)
        feas = feas && all(f -> is_feasible(m, f, n), wp._saf_eq)
        feas = feas && all(f -> is_feasible(m, f, n), wp._sqf_leq)
        feas = feas && all(f -> is_feasible(m, f, n), wp._sqf_eq)
        feas = feas && all(f -> is_feasible(m, f, n), wp._nonlinear_constr)
    end

    if feasible_flag
        interval_objective_used = interval_objective_bound(m, n)
        m._current_xref .= mid(n)
        unsafe_check_fill!(isnan, m._current_xref, 0.0, length(m._current_xref))
    else
        m._lower_objective_value = -Inf
    end
    m._lower_feasibility = feasible_flag
    return
end

"""
$(SIGNATURES)

Runs interval, linear, quadratic contractor methods followed by obbt and a
constraint programming walk up to tolerances specified in
`EAGO.Optimizer` object.
"""
function preprocess!(t::ExtensionType, m::Optimizer)

    reset_relaxation!(m)

    wp = m._working_problem
    params = m._parameters

    # Sets initial feasibility
    feasible_flag = true
    m._obbt_performed_flag = false

    # compute initial volume
    m._initial_volume = prod(upper_variable_bounds(m._current_node) -
                             lower_variable_bounds(m._current_node))
    if params.fbbt_lp_depth >= m._iteration_count
        load_fbbt_buffer!(m)
        for i = 1:m._parameters.fbbt_lp_repetitions
            if feasible_flag
                for j = 1:wp._saf_leq_count
                    !feasible_flag && break
                    saf_leq =  wp._saf_leq[j]
                    feasible_flag &= fbbt!(m, saf_leq)
                end
                !feasible_flag && break

                for j = 1:wp._saf_eq_count
                    !feasible_flag && break
                    saf_eq = wp._saf_eq[j]
                    feasible_flag &= fbbt!(m, saf_eq)
                end
                !feasible_flag && break
            end
        end
        unpack_fbbt_buffer!(m)
    end
    # done after cp to prevent using cp specific flags in cut generation
    set_first_relax_point!(m)

    cp_walk_count = 0
    perform_cp_walk_flag = feasible_flag
    perform_cp_walk_flag &= (params.cp_depth >= m._iteration_count)
    perform_cp_walk_flag &= (cp_walk_count < m._parameters.cp_repetitions)
    while perform_cp_walk_flag
        feasible_flag &= set_constraint_propagation_fbbt!(m)
        !feasible_flag && break
        cp_walk_count += 1
        perform_cp_walk_flag = (cp_walk_count < m._parameters.cp_repetitions)
    end

    obbt_count = 0
    perform_obbt_flag = feasible_flag
    perform_obbt_flag &= (params.obbt_depth >= m._iteration_count)
    perform_obbt_flag &= (obbt_count < m._parameters.obbt_repetitions)

    while perform_obbt_flag
        feasible_flag &= obbt!(m)
        m._obbt_performed_flag = true
        !feasible_flag && break
        obbt_count += 1
        perform_obbt_flag     = (obbt_count < m._parameters.obbt_repetitions)
    end

    m._final_volume = prod(upper_variable_bounds(m._current_node) -
                           lower_variable_bounds(m._current_node))

    m._preprocess_feasibility = feasible_flag

    return
end
preprocess!(m::Optimizer) = preprocess!(m.ext_type, m)

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
    _set_has_value!(wp._objective, false)
    wp._relaxed_evaluator.interval_intersect = false

    if !m._obbt_performed_flag
        relax_constraints!(m, 1)
    end
    relax_objective!(m, 1)

    # Optimizes the object
    d = m.relaxed_optimizer
    MOI.set(d, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    MOI.optimize!(d)

    m._lower_termination_status = MOI.get(d, MOI.TerminationStatus())
    m._lower_primal_status = MOI.get(d, MOI.PrimalStatus())
    m._lower_dual_status = MOI.get(d, MOI.DualStatus())
    status = relaxed_problem_status(m._lower_termination_status,
                                    m._lower_primal_status,
                                    m._lower_dual_status)

    if status == RRS_OPTIMAL
        set_dual!(m)
        m._cut_add_flag = true
        m._lower_feasibility = true
        m._lower_objective_value = MOI.get(d, MOI.ObjectiveValue())
        for i = 1:m._working_problem._variable_count
             m._lower_solution[i] = MOI.get(d, MOI.VariablePrimal(), m._relaxed_variable_index[i])
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

#=
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
        objective_lo = lower_interval_bound(m, m._working_problem._objective, n)
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
=#
