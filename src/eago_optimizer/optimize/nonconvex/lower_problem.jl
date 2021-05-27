
relax_objective!(m, f::Union{SV,SAF}, check_safe) = nothing

function relax_objective!(m, f::BufferedQuadraticIneq, check_safe)
    wp = m._working_problem
    finite_cut_generated = affine_relax_quadratic!(f.func, f.buffer, f.saf,
                                                   m._current_node, m._sol_to_branch_map,
                                                   m._current_xref)
    if finite_cut_generated && (!check_safe || is_safe_cut!(m, f.saf))
        copyto!(wp._objective_saf.terms, f.saf.terms)
        wp._objective_saf.constant = f.saf.constant
        MOI.add(m.relaxed_optimizer, wp._objective_saf, LT_ZERO)
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

    return
end

function reset_relaxation!(m::Optimizer)

    m._cut_iterations = 1
    m._working_problem._relaxed_evaluator.is_first_eval = true
    fill!(m._working_problem._relaxed_evaluator.subexpressions_eval, false)

    m._new_eval_objective = true
    m._new_eval_constraint = true

    delete_nl_constraints!(m)
    delete_objective_cuts!(m)

    return
end

"""
$(FUNCTIONNAME)

"""
function set_first_relax_point!(m::Optimizer)
    m._working_problem._relaxed_evaluator.is_first_eval = true
    m._new_eval_constraint = true
    m._new_eval_objective = true
    m._current_xref .= mid(m._current_node)
    unsafe_check_fill!(isnan, m._current_xref, 0.0, length(m._current_xref))
    return
end

"""
$(TYPEDSIGNATURES)

A routine that adds relaxations for all nonlinear constraints and quadratic constraints
corresponding to the current node to the relaxed problem. This adds an objective cut
(if specified by `objective_cut_on`) and then sets the `_new_eval_constraint` flag
to false indicating that an initial evaluation of the constraints has occurred. If
the `objective_cut_on` flag is `true` then the `_new_eval_objective` flag is also
set to `false` indicating that the objective expression was evaluated.
"""
function relax_all_constraints!(t::ExtensionType, m::Optimizer, k::Int)
    check_safe = (k == 1) ? false : m._parameters.cut_safe_on
    wp = m._working_problem
    wp._relaxed_evaluator.is_first_eval = m._new_eval_constraint
    foreach((i,f) -> relax!(m, f, i, check_safe), enumerate(wp._sqf_leq))
    foreach((i,f) -> relax!(m, f, i, check_safe), enumerate(wp._sqf_eq))
    foreach((i,f) -> relax!(m, f, i, check_safe), enumerate(wp._nonlinear_constr))
    m._new_eval_constraint = false
    objective_cut!(m, check_safe)
    return
end
relax_constraints!(t::ExtensionType, m::Optimizer, q::Int) = relax_all_constraints!(t, m, q)
relax_constraints!(m::Optimizer, k::Int) = relax_constraints!(m.ext_type, m, q)

function relax_problem!(m::Optimizer)
    if m._cut_iterations == 1
        reset_relaxation!(m)
        if m._nonlinear_evaluator_created
            set_node!(wp._relaxed_evaluator, n)
            set_reference_point!(m)
            fill!(wp._relaxed_evaluator.subexpressions_eval, false)
        end
        wp._relaxed_evaluator.is_post = m._parameters.subgrad_tighten
        _set_has_value!(wp._objective, false)
        wp._relaxed_evaluator.interval_intersect = false
        update_relaxed_problem_box!(m)
        set_first_relax_point!(m)
    end
    relax_constraints!(m, k)
    MOI.set(d, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    return
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

function interval_objective_bound(m::Optimizer)
    fL = bound_objective(m)
    if fL > m._lower_objective_value
        m._lower_objective_value = fL
        fill!(m._lower_lvd, 0.0)
        fill!(m._lower_uvd, 0.0)
        m._cut_add_flag = false
    end
    return fL > m._lower_objective_value
end

"""
$(SIGNATURES)

A fallback lower bounding problem that consists of an natural interval extension
calculation. This is called when the optimizer used to compute the lower bound
does not return a termination and primal status code indicating that it
successfully solved the relaxation to a globally optimal point.
"""
function fallback_interval_lower_bound!(m::Optimizer, n::NodeBB)
    feas = true
    wp = m._working_problem
    if !cp_condition(m)
        feas = feas && all(f -> is_feasible(m, f, n), wp._saf_leq)
        feas = feas && all(f -> is_feasible(m, f, n), wp._saf_eq)
        feas = feas && all(f -> is_feasible(m, f, n), wp._sqf_leq)
        feas = feas && all(f -> is_feasible(m, f, n), wp._sqf_eq)
        feas = feas && all(f -> is_feasible(m, f, n), wp._nonlinear_constr)
    end
    if feas
        interval_objective_bound!(m)
        m._current_xref .= mid(n)
        unsafe_check_fill!(isnan, m._current_xref, 0.0, length(m._current_xref))
    else
        m._lower_objective_value = typemin(Float64)
    end
    m._lower_feasibility = feas
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

Checks if a cut should be added and computes a new reference point to add the
cut at. If no cut should be added the constraints not modified in place are
deleted from the relaxed optimizer and the solution is compared with the
interval lower bound. The best lower bound is then used.
"""
function cut_condition(t::ExtensionType, m::Optimizer)
    f_old = m._last_cut_objective
    f_new = m._lower_objective_value
    ϵ_abs = _cut_ϵ_abs(m)
    ϵ_rel = _cut_ϵ_rel(m)
    add_cut_flag = (m._cut_iterations < m._parameters.cut_max_iterations)
    add_cut_flag &= f_new - f_old <= min(ϵ_rel*abs(f_new), ϵ_abs)
    return
end
cut_condition(m::Optimizer) = cut_condition(m.ext_type, m)

"""
$(SIGNATURES)

Constructs and solves the relaxation using the default EAGO relaxation scheme
and optimizer on node `y`.
"""
function lower_problem!(t::ExtensionType, m::Optimizer)
    d = m.relaxed_optimizer

    # start Kelley cutting plane like algorithm
    m._last_cut_objective = typemin(Float64)
    m._lower_objective_value = typemin(Float64)
    set_first_relax_point!(m)
    while true
        relax_problem!(m, k, true)
        m._last_cut_objective = m._lower_objective_value
        MOI.optimize!(d)
        t_status = MOI.get(d, MOI.TerminationStatus())
        p_status = MOI.get(d, MOI.PrimalStatus())
        d_status = MOI.get(d, MOI.DualStatus())
        status = relaxed_problem_status(t_status, p_status, d_status)
        if status != RRS_OPTIMAL
            break
        end
        m._lower_objective_value = MOI.get(d, MOI.ObjectiveValue())
        if cut_condition(m)
            m._cut_iterations += 1
        else
            break
        end
    end

    # activate integrality conditions for MIP & solve MIP subproblem

    t_status = MOI.get(d, MOI.TerminationStatus())
    p_status = MOI.get(d, MOI.PrimalStatus())
    d_status = MOI.get(d, MOI.DualStatus())
    m._lower_termination_status = t_status
    m._lower_primal_status = p_status
    m._lower_dual_status = d_status
    status = relaxed_problem_status(t_status, p_status, d_status)

    if status == RRS_INFEASIBLE
        m._lower_feasibility  = false
        m._lower_objective_value = -Inf
        return
    elseif status == RRS_INVALID
        return fallback_interval_lower_bound!(m, n)
    end

    set_dual!(m)
    m._lower_feasibility = true
    for i = 1:m._working_problem._variable_count
         m._lower_solution[i] = MOI.get(d, MOI.VariablePrimal(), m._relaxed_variable_index[i])
    end
    if status == RRS_DUAL_FEASIBLE
        m._lower_objective_value = MOI.get(d, MOI.DualObjectiveValue())
    end
    interval_objective_bound!(m)
    return
end
lower_problem!(m::Optimizer) = lower_problem!(m.ext_type, m)
