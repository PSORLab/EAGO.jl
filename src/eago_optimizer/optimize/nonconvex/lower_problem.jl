# Copyright (c) 2018: Matthew Wilhelm, Robert Gottlieb, Dimitri Alston,
# Matthew Stuber, and the University of Connecticut (UConn).
# This code is licensed under the MIT license (see LICENSE.md for full details).
################################################################################
# EAGO
# A development environment for robust and global optimization.
# https://github.com/PSORLab/EAGO.jl
################################################################################
# src/eago_optimizer/optimize/nonconvex/lower_problem.jl
# Functions which define how relaxations of subproblems are constructed,
# when domain reductions algorithms are run, and the lower bounding (relaxed
# problem solution subroutines).
################################################################################

"""
$(FUNCTIONNAME)

Add linear objective cut constraint to the `m._subsolvers.relaxed_optimizer`.
"""
function objective_cut!(m::GlobalOptimizer, check_safe::Bool)
    f = m._working_problem._objective_saf
    u = m._global_upper_bound
    if  u < Inf
        b = f.constant
        f.constant = 0.0
        if check_safe && is_safe_cut!(m, f)
            s = LT(u - b + _constraint_tol(m))
            m._affine_objective_cut_ci = MOI.add_constraint(_relaxed_optimizer(m), f, s)
        end
        f.constant = b
        m._new_eval_objective = false
    end
    return
end

"""
    RelaxResultStatus

Status code used internally to determine how to interpret the results from the
solution of a relaxed problem.
"""
@enum(RelaxResultStatus, RRS_OPTIMAL, RRS_DUAL_FEASIBLE, RRS_INFEASIBLE, RRS_INVALID)


"""
$(SIGNATURES)

Take an `MOI.TerminationStatusCode` and two `MOI.ResultStatusCode`s (one each for
the primal and dual status) and return a `RelaxResultStatus`. Returns `RRS_OPTIMAL`
if the codes prove that the subproblem solution was solved to global optimality.
Returns `RRS_INFEASIBLE` if the codes prove that the subproblem solution is
infeasible. Returns `RRS_DUAL_FEASIBLE` if subproblem solution is not optimal
and not proven infeasible, but the dual status is `MOI.FEASIBLE_POINT`. Returns
`RRS_INVALID` otherwise.
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

Update the relaxed constraint by setting the constraint set of `v == x*` ,
`xL_i <= x_i`, and `x_i <= xU_i` for each such constraint added to the relaxed
optimizer. Resets integral valued constraints to either `EqualTo` or `Interval` 
constraints.
"""
function update_relaxed_problem_box!(m::GlobalOptimizer{R,S,Q}) where {R,S,Q<:ExtensionType}
    d = _relaxed_optimizer(m)
    for i = 1:_variable_num(BranchVar(), m)
        l = _lower_bound(BranchVar(), m, i)
        u = _upper_bound(BranchVar(), m, i)
        v = VI(_bvi(m, i))
        if l == u
            ci_vi_et = MOI.add_constraint(d, v, ET(l))
            push!(m._relaxed_variable_et, (ci_vi_et,i))
        else
            ci_vi_lt = MOI.add_constraint(d, v, LT(u))
            ci_vi_gt = MOI.add_constraint(d, v, GT(l))
            m._node_to_sv_leq_ci[i] = ci_vi_lt
            m._node_to_sv_geq_ci[i] = ci_vi_gt
            push!(m._relaxed_variable_lt, (ci_vi_lt,i))
            push!(m._relaxed_variable_gt, (ci_vi_gt,i))
        end
    end
    return
end

const SOLUTION_EPS = 0.05
function store_lower_solution!(m::GlobalOptimizer{R,S,Q}, d::T) where {R,S,Q<:ExtensionType,T}
    for i = 1:_variable_num(FullVar(), m)
        l = _lower_bound(FullVar(), m, i)
        u = _upper_bound(FullVar(), m, i)
        ladj = l + SOLUTION_EPS*(u - l)
        uadj = u - SOLUTION_EPS*(u - l)
        x = MOI.get(d, MOI.VariablePrimal(), m._relaxed_variable_index[i])
        (x < ladj) && (x = ladj)
        (x > uadj) && (x = uadj)
        m._lower_solution[i] = x
    end
    return 
end

function reset_relaxation!(m::GlobalOptimizer{R,S,Q}) where {R,S,Q<:ExtensionType}
    d = _relaxed_optimizer(m)
    m._cut_iterations = 1
    m._obbt_performed_flag = false
    m._working_problem._relaxed_evaluator.is_first_eval = true
    fill!(m._working_problem._relaxed_evaluator.subexpressions_eval, false)

    m._new_eval_objective = true
    m._new_eval_constraint = true

    # Delete added affine constraints
    foreach(c -> MOI.delete(d, c), m._affine_relax_ci)
    empty!(m._affine_relax_ci)

    # Delete variable    
    foreach(c -> MOI.delete(d, c[1]), m._relaxed_variable_et)
    foreach(c -> MOI.delete(d, c[1]), m._relaxed_variable_lt)
    foreach(c -> MOI.delete(d, c[1]), m._relaxed_variable_gt)
    foreach(c -> MOI.delete(d, c), m._relaxed_variable_integer)
    empty!(m._relaxed_variable_et)
    empty!(m._relaxed_variable_lt)
    empty!(m._relaxed_variable_gt)
    empty!(m._relaxed_variable_integer)

    # Delete objective cut
    !isnothing(m._affine_objective_cut_ci) && MOI.delete(d, m._affine_objective_cut_ci)
    return
end

"""
$(TYPEDSIGNATURES)

"""
function set_first_relax_point!(m::GlobalOptimizer)
    if m._cut_iterations == 1
        m._working_problem._relaxed_evaluator.is_first_eval = true
        m._new_eval_constraint = true
        m._new_eval_objective = true
        for i = 1:_variable_num(FullVar(), m)
            l = _lower_bound(FullVar(), m, i)
            u = _upper_bound(FullVar(), m, i)
            if isfinite(l) && isfinite(u)
                x = 0.5*(l + u)
            elseif isfinite(l)
                x = max(0.0, l)
            elseif isfinite(u)
                x = min(0.0, u)
            else
                x = 0.0
            end
            _set_lower_solution!(FullVar(), m, x, i)
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
function relax_all_constraints!(t::ExtensionType, m::GlobalOptimizer, k::Int)
    check_safe = (k == 1) ? false : m._parameters.cut_safe_on
    wp = m._working_problem
    wp._relaxed_evaluator.is_first_eval = m._new_eval_constraint
    foreach(f -> relax!(m, f, k, check_safe), wp._sqf_leq)
    foreach(f -> relax!(m, f, k, check_safe), wp._sqf_eq)
    valid_relax_flag = true
    num_feasible_relax_flag = true
    if valid_relax_flag
        for nl in wp._nonlinear_constr
            valid_cut, feas_cut = relax!(m, nl, k, check_safe)
            valid_relax_flag &= valid_cut
            num_feasible_relax_flag &= feas_cut
        end
    end
    m._new_eval_constraint = false
    (k == 1) && objective_cut!(m, check_safe)
    return valid_relax_flag, num_feasible_relax_flag
end
relax_constraints!(t::ExtensionType, m::GlobalOptimizer, k::Int) = relax_all_constraints!(t, m, k)
relax_constraints!(m::GlobalOptimizer{R,S,Q}, k::Int) where {R,S,Q<:ExtensionType} = relax_constraints!(_ext(m), m, k)

function relax_problem!(m::GlobalOptimizer{R,S,Q}) where {R,S,Q<:ExtensionType}
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
    else
        set_reference_point!(m)
    end
    valid_relax_flag, num_feasible_relax_flag = relax_constraints!(m, m._cut_iterations)
    MOI.set(_relaxed_optimizer(m), MOI.ObjectiveSense(), MOI.MIN_SENSE)
    return valid_relax_flag, num_feasible_relax_flag
end

"""
$(SIGNATURES)

Retrieves the lower and upper duals for variable bounds from the
`relaxed_optimizer` and sets the appropriate values in the
`_lower_lvd` and `_lower_uvd` storage fields.
"""
function set_dual!(m::GlobalOptimizer{R,S,Q}) where {R,S,Q<:ExtensionType}
    d = _relaxed_optimizer(m)
    if MOI.get(d, MOI.DualStatus()) == MOI.FEASIBLE_POINT
        for (c, i) in m._relaxed_variable_lt
            m._lower_uvd[i] = MOI.get(d, MOI.ConstraintDual(), c)
        end
        for (c, i) in m._relaxed_variable_gt
            m._lower_lvd[i] = MOI.get(d, MOI.ConstraintDual(), c)
        end
    else
        fill!(m._lower_lvd, 0.0)
        fill!(m._lower_uvd, 0.0)
    end
    return
end

"""
"""
function interval_objective_bound! end
interval_objective_bound!(m::GlobalOptimizer, f::Nothing, is_first_eval) = nothing
function interval_objective_bound!(m::GlobalOptimizer, f::AffineFunctionIneq, is_first_eval)
    m._working_problem._relaxed_evaluator.is_first_eval = is_first_eval
    fL, fU = bound_objective(m)
    if fL > m._lower_objective_value
        m._lower_objective_value = fL
        fill!(m._lower_lvd, 0.0)
        fill!(m._lower_uvd, 0.0)
        m._cut_add_flag = false
    end
end
function interval_objective_bound!(m::GlobalOptimizer, f, is_first_eval)
    m._working_problem._relaxed_evaluator.is_first_eval = is_first_eval
    if is_first_eval
        m._working_problem._relaxed_evaluator.pass_number = 1
    end
    fL, fU = bound_objective(m)
    fv = _is_input_min(m) ? fL : -fU
    if fv > m._lower_objective_value
        m._lower_objective_value = fv
        fill!(m._lower_lvd, 0.0)
        fill!(m._lower_uvd, 0.0)
        m._cut_add_flag = false
    end
    return
end
interval_objective_bound!(m::GlobalOptimizer, is_first_eval = true) = interval_objective_bound!(m, m._working_problem._objective, is_first_eval)

"""
$(TYPEDSIGNATURES)

Runs contractor methods prior to solving lower bounding problem. By default linear and quadratic 
contractor methods followed by interval constraint propagation then optimization-based bound 
tightening for a specified number of iterations while the subproblem at current node `n` has 
not been proven infeasible.
"""
function preprocess!(t::ExtensionType, m::GlobalOptimizer{R,S,Q}) where {R,S,Q<:ExtensionType}

    feasible_flag = true
    reset_relaxation!(m)
    if _fbbt_lp_depth(m) >= _iteration_count(m)
        load_fbbt_buffer!(m)
        for _ = 1:_fbbt_lp_repetitions(m)
            ns = NodeBB(_current_node(m))
            for f in m._working_problem._saf_leq
                !(feasible_flag = feasible_flag && fbbt!(m, f)) && break
            end
            !feasible_flag && break
            for f in m._working_problem._saf_eq
                !(feasible_flag = feasible_flag && fbbt!(m, f)) && break
            end
            (same_box(ns,_current_node(m),0.0) || !feasible_flag) && break
        end
        unpack_fbbt_buffer!(m)
    end

    # Done after CP to prevent using CP specific flags in cut generation
    set_first_relax_point!(m)
    # Nonlinear CP can detect infeasibility and bound objective even if
    # the relaxation is ill-posed, so one is always used to mitigate numerical issues 
    cp_reps = _cp_depth(m) >= _iteration_count(m) ? _cp_repetitions(m) : 0
    for _ = 1:cp_reps
        ns = NodeBB(_current_node(m))
        feasible_flag = feasible_flag && set_constraint_propagation_fbbt!(m)
        (same_box(ns,_current_node(m),0.0) || !feasible_flag) && break
    end

    if _obbt_depth(m) >= _iteration_count(m)
        for k = 1:_obbt_repetitions(m) 
            ns = NodeBB(_current_node(m))
            feasible_flag = feasible_flag && obbt!(m)
            m._obbt_performed_flag = true
            (same_box(ns,_current_node(m),0.0) || !feasible_flag) && break
        end
    end

    m._preprocess_feasibility = feasible_flag
    return
end
preprocess!(m::GlobalOptimizer{R,S,Q}) where {R,S,Q<:ExtensionType} = preprocess!(_ext(m), m)

"""
$(TYPEDSIGNATURES)

Returns `true` if a cut should be added and computes a new reference point to add the
cut at. By default, checks that `cut_max_iterations` are not exceeded and that the 
improvement in the objective value associated with the previous cut is greater than
both an absolute tolerance `cut_ϵ_abs` and a relative tolerance `cut_ϵ_rel`. Returns
`false` otherwise.
"""
function cut_condition(t::ExtensionType, m::GlobalOptimizer)
    obj_old = m._last_cut_objective
    obj_new = m._lower_objective_value
    flag = m._cut_iterations < _cut_max_iterations(m)
    flag &= obj_new - obj_old > _cut_ϵ_abs(m)
    flag &= obj_new - obj_old > _cut_ϵ_rel(m)*abs(obj_new)
    return flag
end
cut_condition(m::GlobalOptimizer{R,S,Q}) where {R,S,Q<:ExtensionType} = cut_condition(_ext(m), m)

"""
$(SIGNATURES)

Returns `true` that the subproblem at the current node `n` has participating integer variables
that have not been fixed to constant valued as the branch-and-bound algorithm progresses. Returns
`false` otherwise.
"""
is_integer_subproblem(m) = !continuous(_current_node(m))

"""
$(TYPEDSIGNATURES)

Constructs a relaxation of the MINLP on node `y` and solves it using the default EAGO 
relaxation scheme. By default, EAGO applies Kelley's algorithm (from Kelley Jr., J.E.: 
The cutting-plane method for solving convex programs. J. Soc. Ind. Appl. Math. 8(4), 
703 to 712 (1960)) while `cut_condition(m)` returns `true` then activates the integrality
constraints of the relaxed problems and solves the resulting MILP relaxation. results
are stored to the `_lower_solution`, `_lower_termination_status`, `_lower_primal_status`,
`_lower_dual_status`, `_lower_objective_value`, and `_lower_feasibility`. Further, lower
and upper variable duals are stored `_lower_lvd` and `_lower_uvd`, respectively, for use
in duality based bound tightening. If relaxation-based bounds are weaker or cutting-planes 
are numerically poorly ill-posed, then interval bounds are used instead. If the problem is
dual feasible but the primal status is ambiguous the dual objective value is used for the 
lower bound to avoid numerical issues.
"""
function lower_problem!(t::ExtensionType, m::GlobalOptimizer{R,S,Q}) where {R,S,Q<:ExtensionType}

    num_feasible_relax_flag = true

    d = _relaxed_optimizer(m)
    m._last_cut_objective = typemin(Float64)
    m._lower_objective_value = typemin(Float64)

    t_status = MOI.OPTIMIZE_NOT_CALLED
    p_status = MOI.OTHER_RESULT_STATUS
    d_status = MOI.OTHER_RESULT_STATUS
    status = RRS_INVALID

    set_first_relax_point!(m)
    MOI.set(d, MOI.ObjectiveFunction{SAF}(), m._working_problem._objective_saf)

    while true
        valid_prob, feas_flag = relax_problem!(m)
        if !feas_flag 
            num_feasible_relax_flag = false
            break
        end
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
            store_lower_solution!(m, d)
            m._cut_iterations += 1
        else
            break
        end
    end
    if !num_feasible_relax_flag
        status = RRS_INFEASIBLE
    end

    # Activate integrality conditions for MIP and solve MIP subproblem
    if is_integer_subproblem(m) && (status !== RRS_INFEASIBLE)
        m._last_cut_objective = m._lower_objective_value
        for i = 1:_variable_num(BranchVar(), m)
            l = _lower_bound(BranchVar(), m, i)
            u = _upper_bound(BranchVar(), m, i)
            if is_integer(BranchVar(), m, i) && (l != u)
                c_integer = MOI.add_constraint(d, VI(_bvi(m, i)), INT())
                push!(m._relaxed_variable_integer, c_integer)
            end
        end
        MOI.optimize!(d)
        t_status = MOI.get(d, MOI.TerminationStatus())
        p_status = MOI.get(d, MOI.PrimalStatus())
        d_status = MOI.get(d, MOI.DualStatus())
        status = relaxed_problem_status(t_status, p_status, d_status)
        if status == RRS_OPTIMAL
            m._lower_objective_value = MOI.get(d, MOI.ObjectiveValue())
        end
    end

    # Check status, if not feasible/infeasible then fallback to interval bounds
    if status == RRS_OPTIMAL
        m._lower_objective_value = MOI.get(d, MOI.ObjectiveValue())
    end

    m._lower_termination_status = t_status
    m._lower_primal_status = p_status
    m._lower_dual_status = d_status
    status = relaxed_problem_status(t_status, p_status, d_status)
    if !num_feasible_relax_flag
        status = RRS_INFEASIBLE
    end
    if status == RRS_INFEASIBLE
        m._lower_feasibility  = false
        m._lower_objective_value = -Inf
        return
    end

    # Set dual values
    set_dual!(m)
    m._lower_feasibility = true
    store_lower_solution!(m, d)
    if status == RRS_DUAL_FEASIBLE
        m._lower_objective_value = MOI.get(d, MOI.DualObjectiveValue())
    end
    interval_objective_bound!(m, true) 
    return
end
lower_problem!(m::GlobalOptimizer{R,S,Q}) where {R,S,Q<:ExtensionType} = lower_problem!(_ext(m), m)