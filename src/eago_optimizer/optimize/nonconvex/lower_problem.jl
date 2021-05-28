"""
$(FUNCTIONNAME)

Adds linear objective cut constraint to the `x.relaxed_optimizer`.
"""
function objective_cut!(m::Optimizer, check_safe::Bool)
    wp = m._working_problem
    f = wp._objective_saf
    u = m._global_upper_bound
    if  u < Inf
        b = f.constant
        f.constant = 0.0
        if check_safe && is_safe_cut!(m, f)
            s = LT(u - b + _constraint_tol(m))
            c = MOI.add_constraint(m.relaxed_optimizer, f, s)::CI{SAF,LT}
            m._affine_objective_cut_ci = c
        end
        f.constant = b
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

    # delete added affine constraints
    foreach(c -> MOI.delete(m.relaxed_optimizer, c)::Nothing, m._affine_relax_ci)
    empty!(m._affine_relax_ci)

    # delete objective cut
    if m._affine_objective_cut_ci !== nothing
        MOI.delete(m.relaxed_optimizer, m._affine_objective_cut_ci)::Nothing
    end

    return
end

val_or_zero(x) = isnan(x) ? 0.0 : x
"""
$(FUNCTIONNAME)

"""
function set_first_relax_point!(m::Optimizer)
    if m._cut_iterations == 1
        m._working_problem._relaxed_evaluator.is_first_eval = true
        m._new_eval_constraint = true
        m._new_eval_objective = true
        for i = 1:_variable_num(BranchVar(), m)
            x = _mid(BranchVar(), m, i)
            _set_lower_solution!(BranchVar(), m, val_or_zero(x), i)
        end
    end
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
    foreach(x -> relax!(m, x, k, check_safe), wp._sqf_leq)
    foreach(x -> relax!(m, x, k, check_safe), wp._sqf_eq)
    foreach(x -> relax!(m, x, k, check_safe),wp._nonlinear_constr)
    relax!(m, wp._objective, k, check_safe)
    m._new_eval_constraint = false
    (k == 1) && objective_cut!(m, check_safe)
    return
end
relax_constraints!(t::ExtensionType, m::Optimizer, k::Int) = relax_all_constraints!(t, m, k)
relax_constraints!(m::Optimizer, k::Int) = relax_constraints!(m.ext_type, m, k)

function relax_problem!(m::Optimizer)
    wp = m._working_problem
    if m._cut_iterations == 1
        reset_relaxation!(m)
        if m._nonlinear_evaluator_created
            set_node!(wp._relaxed_evaluator, m._current_node)
            set_reference_point!(m)
            fill!(wp._relaxed_evaluator.subexpressions_eval, false)
        end
        wp._relaxed_evaluator.is_post = m._parameters.subgrad_tighten
        _set_has_value!(wp._objective, false)
        wp._relaxed_evaluator.interval_intersect = false
        update_relaxed_problem_box!(m)
        set_first_relax_point!(m)
    end
    relax_constraints!(m, m._cut_iterations)
    MOI.set(m.relaxed_optimizer, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    @show "ran relaxed problem"
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

function interval_objective_bound!(m::Optimizer)
    fL = bound_objective(m)*m._obj_mult
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
            for f in wp._saf_leq
                !(feasible_flag = feasible_flag && fbbt!(m, f)) && break
            end
            !feasible_flag && break
            for f in wp._saf_eq
                !(feasible_flag = feasible_flag && fbbt!(m, f)) && break
            end
            !feasible_flag && break
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
        perform_obbt_flag = (obbt_count < m._parameters.obbt_repetitions)
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
    obj_old = m._last_cut_objective
    obj_new = m._lower_objective_value
    flag = m._cut_iterations < _cut_max_iterations(m)
    flag &= obj_new - obj_old > _cut_ϵ_abs(m)
    flag &= obj_new - obj_old > _cut_ϵ_rel(m)*abs(obj_new)
    return flag
end
cut_condition(m::Optimizer) = cut_condition(m.ext_type, m)

"""
$(SIGNATURES)

Constructs and solves the relaxation using the default EAGO relaxation scheme
and optimizer on node `y`.
"""
function lower_problem!(t::ExtensionType, m::Optimizer)
    m._last_cut_objective = typemin(Float64)
    m._lower_objective_value = typemin(Float64)
    set_first_relax_point!(m)
    MOI.set(m.relaxed_optimizer, MOI.ObjectiveFunction{SAF}(), m._working_problem._objective_saf)
    while true
        relax_problem!(m)
        m._last_cut_objective = m._lower_objective_value
        MOI.optimize!(m.relaxed_optimizer)::Nothing
        t_status = MOI.get(m.relaxed_optimizer, MOI.TerminationStatus())::MOI.TerminationStatusCode
        p_status = MOI.get(m.relaxed_optimizer, MOI.PrimalStatus())::MOI.ResultStatusCode
        d_status = MOI.get(m.relaxed_optimizer, MOI.DualStatus())::MOI.ResultStatusCode
        status = relaxed_problem_status(t_status, p_status, d_status)
        if status != RRS_OPTIMAL
            break
        end
        m._lower_objective_value = MOI.get(m.relaxed_optimizer, MOI.ObjectiveValue())::Float64
        m._lower_objective_value *= m._obj_mult
        @show m._cut_iterations
        @show m._lower_objective_value
        if cut_condition(m)
            m._cut_iterations += 1
        else
            break
        end
    end

    # activate integrality conditions for MIP & solve MIP subproblem

    t_status = MOI.get(m.relaxed_optimizer, MOI.TerminationStatus())::MOI.TerminationStatusCode
    p_status = MOI.get(m.relaxed_optimizer, MOI.PrimalStatus())::MOI.ResultStatusCode
    d_status = MOI.get(m.relaxed_optimizer, MOI.DualStatus())::MOI.ResultStatusCode
    m._lower_termination_status = t_status
    m._lower_primal_status = p_status
    m._lower_dual_status = d_status
    status = relaxed_problem_status(t_status, p_status, d_status)
    @show MOI.get(m.relaxed_optimizer, MOI.ObjectiveValue())
    @show t_status
    @show p_status
    @show d_status

    if status == RRS_INFEASIBLE
        m._lower_feasibility  = false
        m._lower_objective_value = -Inf
        return
    elseif status == RRS_INVALID
        return fallback_interval_lower_bound!(m, n)::Nothing
    end

    set_dual!(m)
    m._lower_feasibility = true
    for i = 1:m._working_problem._variable_count
         m._lower_solution[i] = MOI.get(m.relaxed_optimizer, MOI.VariablePrimal(), m._relaxed_variable_index[i])::Float64
    end
    if status == RRS_DUAL_FEASIBLE
        m._lower_objective_value = MOI.get(m.relaxed_optimizer, MOI.DualObjectiveValue())::Float64
        m._lower_objective_value *= m._obj_mult
    end
    interval_objective_bound!(m)
    return
end
lower_problem!(m::Optimizer) = lower_problem!(m.ext_type, m)
