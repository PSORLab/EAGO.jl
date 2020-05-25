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
# src/eago_optimizer/subroutines.jl
# Default subroutines for EAGO's global optimizer.
#############################################################################

"""
$(SIGNATURES)

Selects node with the lowest lower bound in stack.
"""
function node_selection!(t::ExtensionType, x::Optimizer)
    x._node_count -= 1
    x._current_node = popmin!(x._stack)
    return
end

"""
$(SIGNATURES)

Creates two nodes from current_node using information available the `x`
and stores them to the stack. By default, relative width bisection is perfomed
at a point `branch_pnt` which is a convex combination
(parameter: branch_cvx_factor) of the solution to the relaxation and
the midpoint of the node. If this solution lies within `branch_offset/width` of
a bound then the branch point is moved to a distance of `branch_offset/width`
from the bound.
"""
function branch_node!(t::ExtensionType, m::Optimizer)

    y = m._current_node
    lvbs = y.lower_variable_bounds
    uvbs = y.upper_variable_bounds

    max_pos = 0
    max_val = -Inf
    temp_max = 0.0

    flag = true
    for i = 1:m._variable_number
        flag = @inbounds ~m._fixed_variable[i]
        flag &= @inbounds m.branch_variable[i]
        if flag
            vi = @inbounds m._variable_info[i]
            temp_max = @inbounds uvbs[i] - lvbs[i]
            temp_max /= vi.upper_bound - vi.lower_bound
            if temp_max > max_val
                max_pos = i
                max_val = temp_max
            end
        end
    end

    lvb = @inbounds lvbs[max_pos]
    uvb = @inbounds uvbs[max_pos]
    lsol = @inbounds m._lower_solution[max_pos]
    cvx_f = m._parameters.branch_cvx_factor
    cvx_g = m._parameters.branch_offset
    branch_pnt = cvx_f*lsol + (1.0 - cvx_f)*(lvb + uvb)/2.0
    if branch_pnt < lvb*(1.0 - cvx_g) + cvx_g*uvb
        branch_pnt = (1.0 - cvx_g)*lvb + cvx_g*uvb
    elseif branch_pnt > cvx_g*lvb + (1.0 - cvx_g)*uvb
        branch_pnt = cvx_g*lvb + (1.0 - cvx_g)*uvb
    end
    N1::Interval{Float64} = Interval{Float64}(lvb, branch_pnt)
    N2::Interval{Float64} = Interval{Float64}(branch_pnt, uvb)
    lvb_1 = copy(lvbs)
    uvb_1 = copy(uvbs)
    lvb_2 = copy(lvbs)
    uvb_2 = copy(uvbs)
    @inbounds lvb_1[max_pos] = N1.lo
    @inbounds uvb_1[max_pos] = N1.hi
    @inbounds lvb_2[max_pos] = N2.lo
    @inbounds uvb_2[max_pos] = N2.hi

    lower_bound = max(y.lower_bound, m._lower_objective_value)
    upper_bound = min(y.upper_bound, m._upper_objective_value)
    m._maximum_node_id += 1
    X1 = NodeBB(lvb_1, uvb_1, lower_bound, upper_bound, y.depth + 1, m._maximum_node_id)
    m._maximum_node_id += 1
    X2 = NodeBB(lvb_2, uvb_2, lower_bound, upper_bound, y.depth + 1, m._maximum_node_id)

    push!(m._stack, X1)
    push!(m._stack, X2)

    m._node_repetitions = 1
    m._node_count += 2

    return
end

"""
$(SIGNATURES)

Stores the current node to the stack after updating lower/upper bounds.
"""
function single_storage!(t::ExtensionType, m::Optimizer)
    y = m._current_node
    m._node_repetitions += 1
    m._node_count += 1
    lower_bound = max(y.lower_bound, m._lower_objective_value)
    upper_bound = min(y.upper_bound, m._upper_objective_value)
    push!(m._stack, NodeBB(y.lower_variable_bounds, y.upper_variable_bounds,
                           lower_bound, upper_bound, y.depth, y.id))
    return
end

"""
$(SIGNATURES)

Selects and deletes nodes from stack with lower bounds greater than global
upper bound.
"""
function fathom!(t::ExtensionType, m::Optimizer)
    upper = m._global_upper_bound
    continue_flag = ~isempty(m._stack)
    while continue_flag
        max_node = maximum(m._stack)
        max_check = (max_node.lower_bound > upper)
        if max_check
            popmax!(m._stack)
            m._node_count -= 1
            if isempty(m._stack)
                continue_flag = false
            end
        else
            if ~max_check
                continue_flag = false
            elseif isempty(m._stack)
                continue_flag = false
            end
        end
    end
    return
end

"""
$(SIGNATURES)

Checks to see if current node should be reprocessed.
"""
repeat_check(t::ExtensionType, m::Optimizer) = false

relative_gap(L::Float64, U::Float64) = ((L > -Inf) && (U < Inf)) ?  abs(U - L)/(max(abs(L), abs(U))) : Inf
relative_tolerance(L::Float64, U::Float64, tol::Float64) = relative_gap(L, U)  > tol || ~(L > -Inf)

"""
$(SIGNATURES)

Checks for termination of algorithm due to satisfying absolute or relative
tolerance, infeasibility, or a specified limit, returns a boolean valued true
if algorithm should continue.
"""
function termination_check(t::ExtensionType, m::Optimizer)

    node_in_stack = length(m._stack)
    L = m._global_lower_bound
    U = m._global_upper_bound

    if node_in_stack === 0

        if m._first_solution_node > 0
            m._termination_status_code = MOI.OPTIMAL
            m._result_status_code = MOI.FEASIBLE_POINT
            (m._parameters.verbosity >= 3) && println("Empty Stack: Exhaustive Search Finished")
        else
            m._termination_status_code = MOI.INFEASIBLE
            m._result_status_code = MOI.INFEASIBILITY_CERTIFICATE
            (m._parameters.verbosity >= 3) && println("Empty Stack: Infeasible")
        end

    elseif node_in_stack >= m._parameters.node_limit

        m._termination_status_code = MOI.NODE_LIMIT
        m._result_status_code = MOI.UNKNOWN_RESULT_STATUS
        (m._parameters.verbosity >= 3) && println("Node Limit Exceeded")

    elseif m._iteration_count >= m._parameters.iteration_limit

        m._termination_status_code = MOI.ITERATION_LIMIT
        m._result_status_code = MOI.UNKNOWN_RESULT_STATUS
        (m._parameters.verbosity >= 3) && println("Maximum Iteration Exceeded")

    elseif ~relative_tolerance(L, U, m._parameters.relative_tolerance)

        m._termination_status_code = MOI.OPTIMAL
        m._result_status_code = MOI.FEASIBLE_POINT
        (m._parameters.verbosity >= 3) && println("Relative Tolerance Achieved")

    elseif (U - L) < m._parameters.absolute_tolerance

        m._termination_status_code = MOI.OPTIMAL
        m._result_status_code = MOI.FEASIBLE_POINT
        (m._parameters.verbosity >= 3) && println("Absolute Tolerance Achieved")

    elseif m._run_time > m._parameters.time_limit

        m._termination_status_code = MOI.TIME_LIMIT
        m._result_status_code = MOI.UNKNOWN_RESULT_STATUS
        (m._parameters.verbosity >= 3) && println("Time Limit Exceeded")

    else

        return false

    end

    return true
end

"""
$(SIGNATURES)

Checks for convergence of algorithm with respect to absolute and/or relative
tolerances.
"""
function convergence_check(t::ExtensionType, m::Optimizer)

  L = m._lower_objective_value
  U = m._global_upper_bound
  t = (U - L) <= m._parameters.absolute_tolerance
  if (U < Inf) & (L > Inf)
      t |= (abs(U - L)/(max(abs(L), abs(U))) <= m._parameters.relative_tolerance)
  end
  return t
end

"""
$(SIGNATURES)

Takes an `MOI.TerminationStatusCode` and a `MOI.ResultStatusCode` and returns
the tuple `(valid_result::Bool, feasible::Bool)`. The value `valid_result` is
`true` if the pair of codes prove that either the subproblem solution was solved
to global optimality or the subproblem solution is infeasible. The value of
`feasible` is true if the problem is feasible and false if the problem is infeasible.
"""
function is_globally_optimal(t::MOI.TerminationStatusCode, r::MOI.ResultStatusCode)

    feasible = false
    valid_result = false

    if (t === MOI.INFEASIBLE && r == MOI.INFEASIBILITY_CERTIFICATE)
        valid_result = true
    elseif (t === MOI.INFEASIBLE && r === MOI.NO_SOLUTION)
        valid_result = true
    elseif (t === MOI.INFEASIBLE && r === MOI.UNKNOWN_RESULT_STATUS)
        valid_result = true
    elseif (t === MOI.OPTIMAL && r === MOI.FEASIBLE_POINT)
        valid_result = true
        feasible = true
    elseif (t === MOI.INFEASIBLE_OR_UNBOUNDED && r === MOI.NO_SOLUTION)
        valid_result = true
        feasible = false
    end

    return valid_result, feasible
end

"""
$(SIGNATURES)

Takes an `MOI.TerminationStatusCode` and a `MOI.ResultStatusCode` and returns `true`
if this corresponds to a solution that is proven to be feasible.
Returns `false` otherwise.
"""
function is_feasible_solution(t::MOI.TerminationStatusCode, r::MOI.ResultStatusCode)

    termination_flag = false
    result_flag = false

    (t === MOI.OPTIMAL) && (termination_flag = true)
    (t === MOI.LOCALLY_SOLVED) && (termination_flag = true)

    # This is default solver specific... the acceptable constraint tolerances
    # are set to the same values as the basic tolerance. As a result, an
    # acceptably solved solution is feasible but non necessarily optimal
    # so it should be treated as a feasible point
    if (t === MOI.ALMOST_LOCALLY_SOLVED) && (r === MOI.NEARLY_FEASIBLE_POINT)
        termination_flag = true
        result_flag = true
    end

    (r === MOI.FEASIBLE_POINT) && (result_flag = true)

    return (termination_flag && result_flag)
end

"""
$(SIGNATURES)

Retrieves the lower and upper duals for variable bounds from the
`relaxed_optimizer` and sets the appropriate values in the
`_lower_lvd` and `_lower_uvd` storage fields.
"""
function set_dual!(m::Optimizer)

    opt = m.relaxed_optimizer
    lower_variable_lt_indx = m._lower_variable_lt_indx
    lower_variable_lt = m._lower_variable_lt
    lower_variable_gt_indx = m._lower_variable_gt_indx
    lower_variable_gt = m._lower_variable_gt

    for i = 1:length(lower_variable_lt_indx)
        @inbounds m._lower_uvd[@inbounds lower_variable_lt_indx[i]] = MOI.get(opt, MOI.ConstraintDual(), @inbounds lower_variable_lt[i])
    end
    for i = 1:length(lower_variable_gt_indx)
        @inbounds m._lower_lvd[@inbounds lower_variable_gt_indx[i]] = MOI.get(opt, MOI.ConstraintDual(), @inbounds lower_variable_gt[i])
    end

    return
end

"""
$(SIGNATURES)

Runs interval, linear, quadratic contractor methods followed by obbt and a
constraint programming walk up to tolerances specified in
`EAGO.Optimizer` object.
"""
function preprocess!(t::ExtensionType, m::Optimizer)

    # Sets initial feasibility
    feasible_flag = true
    rept = 0

    m._initial_volume = prod(upper_variable_bounds(m._current_node) -
                             lower_variable_bounds(m._current_node))

    # runs poor man's LP contractor
    if (m._parameters.lp_depth >= m._iteration_count) && feasible_flag
        for func in
        end
    end

    m._obbt_performed_flag = false
    if (m._parameters.obbt_depth >= m._iteration_count) && feasible_flag
        m._obbt_performed_flag = true
        for i = 1:m._parameters.obbt_repetitions
            feasible_flag = obbt(m)
            !feasible_flag && break
        end
    end

    if (m._parameters.cp_depth >= m._iteration_count) && feasible_flag
        feasible_flag = cpwalk(m)
    end

    m._final_volume = prod(upper_variable_bounds(m._current_node) -
                           lower_variable_bounds(m._current_node))
    m._preprocess_feasibility = feasible_flag

    return
end

"""
$(SIGNATURES)

Updates the relaxed constraint by setting the constraint set of `v == x*`` ,
`xL_i <= x_i`, and `x_i <= xU_i` for each such constraint added to the relaxed
optimizer.
"""
function update_relaxed_problem_box!(m::Optimizer, y::NodeBB)

    opt = m.relaxed_optimizer
    lower_node_bnd = y.lower_variable_bounds
    upper_node_bnd = y.upper_variable_bounds

    lower_variable_et = m._lower_variable_et
    lower_variable_et_indx = m._lower_variable_et_indx
    for i = 1:length(lower_variable_et_indx)
        @inbounds ci = lower_variable_et[i]
        @inbounds ni = lower_variable_et_indx[i]
        @inbounds vb = lower_node_bnd[ni]
        MOI.set(opt, MOI.ConstraintSet(), ci, ET(vb))
    end

    lower_variable_lt = m._lower_variable_lt
    lower_variable_lt_indx = m._lower_variable_lt_indx
    for i = 1:length(lower_variable_lt_indx)
        @inbounds ci = lower_variable_lt[i]
        @inbounds ni = lower_variable_lt_indx[i]
        @inbounds vb = upper_node_bnd[ni]
        MOI.set(opt, MOI.ConstraintSet(), ci, LT(vb))
    end

    lower_variable_gt = m._lower_variable_gt
    lower_variable_gt_indx = m._lower_variable_gt_indx
    for i = 1:length(x._lower_variable_gt_indx)
        @inbounds ci = lower_variable_gt[i]
        @inbounds ni = lower_variable_gt_indx[i]
        @inbounds vb = lower_node_bnd[ni]
        MOI.set(opt, MOI.ConstraintSet(), ci, GT(vb))
    end

    return
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

    if !cp_condition(m)
        for i = 1:m._working_problem._saf_leq_count
            saf_leq = @inbounds m._working_problem._saf_leq[i]
            feasible_flag &= (lower_interval_bound(saf_leq, n) <= 0.0)
            !feasible_flag && break
        end

        if feasible_flag
            for i = 1:m._working_problem._saf_eq_count
                saf_eq = @inbounds m._working_problem._saf_eq[i]
                lower_value, upper_value = interval_bound(saf_eq, n)
                feasible_flag &= (lower_value <= 0.0 <= upper_value)
                !feasible_flag && break
            end
        end

        if feasible_flag
            for i = 1:m._working_problem._sqf_leq_count
                sqf_leq = @inbounds m._working_problem._saf_leq[i]
                feasible_flag &= (lower_interval_bound(sqf_leq, n) <= 0.0)
                !feasible_flag && break
            end
        end

        if feasible_flag
            for i = 1:m._working_problem._sqf_eq_count
                sqf_eq = @inbounds m._working_problem._sqf_eq[i]
                lower_value, upper_value = interval_bound(sqf_eq, n)
                feasible_flag &= (lower_value <= 0.0 <= upper_value)
                !feasible_flag && break
            end
        end
    end

    if feasible_flag
        interval_objective_used = interval_objective_bound(m, n)
    else
        x._lower_objective_value = -Inf
    end
    x._lower_feasibility = feasible_flag

    return
end

"""
$(SIGNATURES)

Constructs and solves the relaxation using the default EAGO relaxation scheme
and optimizer on node `y`.
"""
function lower_problem!(t::ExtensionType, m::Optimizer)

    n = m._current_node

    if ~m._obbt_performed_flag
        @. m._current_xref = 0.5*(n.lower_variable_bounds + n.upper_variable_bounds)
        unsafe_check_fill!(isnan, m._current_xref, 0.0, length(m._current_xref))
        update_relaxed_problem_box!(m, n)
        relax_problem!(m, m._current_xref, 1)
    end

    relax_objective!(m, m._current_xref)
    if m.objective_cut_on
        objective_cut_linear!(m, 1)
    end

    # Optimizes the object
    opt = m.relaxed_optimizer
    MOI.optimize!(opt)

    m._lower_termination_status = MOI.get(opt, MOI.TerminationStatus())
    m._lower_result_status = MOI.get(opt, MOI.PrimalStatus())
    valid_flag, feas_flag = is_globally_optimal(m._lower_termination_status, m._lower_result_status)

    if valid_flag
        if feas_flag
            m._lower_feasibility = true
            m._lower_objective_value = MOI.get(opt, MOI.ObjectiveValue())
            @inbounds m._lower_solution[:] = MOI.get(opt, MOI.VariablePrimal(), m._lower_variable_index)
            interval_objective_used = interval_objective_bound(m, n)
            if !interval_objective_used
                set_dual!(m)
                m._cut_add_flag = true
            end
        else
            m._cut_add_flag = false
            m._lower_feasibility  = false
            m._lower_objective_value = -Inf
        end
    else
        fallback_interval_lower_bound!(m, n)
    end
    return
end

"""
$(SIGNATURES)

Updates the internal storage in the optimizer after a valid feasible cut is added.
"""
function cut_update(x::Optimizer)

    x._cut_feasibility = true

    opt = x.relaxed_optimizer
    obj_val = MOI.get(opt, MOI.ObjectiveValue())
    prior_obj_val = (x._cut_iterations == 2) ? x._lower_objective_value : x._cut_objective_value

    if prior_obj_val < obj_val
        x._cut_objective_value = obj_val
        x._lower_objective_value = obj_val
        x._lower_termination_status = x._cut_termination_status
        x._lower_result_status = x._cut_result_status
        @inbounds x._cut_solution[:] = MOI.get(opt, MOI.VariablePrimal(), x._lower_variable_index)
        copyto!(x._lower_solution, x._cut_solution)
        set_dual!(x)
        x._cut_add_flag = true
    else
        x._cut_add_flag = false
    end

    return
end


"""
$(SIGNATURES)

Checks if a cut should be added and computes a new reference point to add the
cut at. If no cut should be added the constraints not modified in place are
deleted from the relaxed optimizer and the solution is compared with the
interval lower bound. The best lower bound is then used.
"""
function cut_condition(t::ExtensionType, x::Optimizer)

    flag = x._cut_add_flag
    flag &= (x._cut_iterations < x.cut_max_iterations)
    y = x._current_node

    if flag
        xprior = x._current_xref
        xsol = (x._cut_iterations > 1) ? x._cut_solution : x._lower_solution
        xnew = (1.0 - x.cut_cvx)*mid(y) + x.cut_cvx*xsol
        if norm((xprior - xnew)./diam(y), 1) > x.cut_tolerance
            x._current_xref = xnew
            flag &= true
        else
            flag &= false
        end
    end

    if ~flag
        # if not further cuts then empty the added cuts
        for i in 2:x._cut_iterations
            for ci in x._quadratic_ci_leq[i]
                MOI.delete(x.relaxed_optimizer, ci)
            end
            for ci in x._quadratic_ci_geq[i]
                MOI.delete(x.relaxed_optimizer, ci)
            end
            for (ci1,ci2) in x._quadratic_ci_eq[i]
                MOI.delete(x.relaxed_optimizer, ci1)
                MOI.delete(x.relaxed_optimizer, ci2)
            end
            for ci in x._lower_nlp_affine[i]
                MOI.delete(x.relaxed_optimizer, ci)
            end
            for ci in x._upper_nlp_affine[i]
                MOI.delete(x.relaxed_optimizer, ci)
            end
        end
        if x._objective_cut_set !== -1
            if x._objective_type === SINGLE_VARIABLE
                obj_indx = x._objective_sv.variable.value
                objective_lo = @inbounds y.lower_variable_bounds[obj_indx]
            elseif x._objective_type === SCALAR_AFFINE
                objective_lo = interval_bound(x._objective_saf, y, true)
            elseif x._objective_type === SCALAR_QUADRATIC
                    objective_lo = interval_bound(x._objective_sqf, y, true)
            end
        end
    end

    # check to see if interval bound is preferable
    if x._lower_feasibility
        if x._objective_type === NONLINEAR
            objective_lo = eval_objective_lo(x._relaxed_evaluator)
        elseif x._objective_type === SINGLE_VARIABLE
                obj_indx = x._objective_sv.variable.value
                objective_lo = @inbounds y.lower_variable_bounds[obj_indx]
        elseif x._objective_type === SCALAR_AFFINE
                objective_lo = interval_bound(x._objective_saf, y, true)
        elseif x._objective_type === SCALAR_QUADRATIC
                objective_lo = interval_bound(x._objective_sqf, y, true)
        end
        if objective_lo > x._lower_objective_value
            x._lower_objective_value = objective_lo
            fill!(x._lower_lvd, 0.0)
            fill!(x._lower_uvd, 0.0)
        end
    end
    x._cut_iterations += 1

    return flag
end

"""
$(SIGNATURES)

Adds a cut for each constraint and the objective function to the subproblem.
"""
function add_cut!(t::ExtensionType, x::Optimizer)

    relax_problem!(x, x._current_xref, x._cut_iterations)
    relax_objective!(x, x._current_xref)
    if x.objective_cut_on
        objective_cut_linear!(x, x._cut_iterations)
    end

    # Optimizes the object
    opt = x.relaxed_optimizer
    MOI.optimize!(opt)

    x._cut_termination_status = MOI.get(opt, MOI.TerminationStatus())
    x._cut_result_status = MOI.get(opt, MOI.PrimalStatus())
    valid_flag, feas_flag = is_globally_optimal(x._cut_termination_status, x._cut_result_status)

    if valid_flag
        if feas_flag
            cut_update(x)
        else
            x._cut_add_flag = false
            x._lower_feasibility  = false
            x._lower_objective_value = -Inf
        end
    else
        x._cut_add_flag = false
    end

    return
end

"""
$(SIGNATURES)

Default check to see if the upper bounding problem should be run. By default,
The upper bounding problem is run on every node up to depth `upper_bounding_depth`
and is triggered with a probability of `0.5^(depth - upper_bounding_depth)`
afterwards.
"""
function default_nlp_heurestic(x::Optimizer, y::NodeBB)
    bool = false
    ubd_limit = x._parameters.upper_bounding_depth
    depth = y.depth
    bool |= (depth <= ubd_limit)
    bool |= (rand() < 0.5^(depth - x._parameters.upper_bounding_depth))
    return bool
end

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
function solve_local_nlp!(x::Optimizer)

    y = x._current_node

    if default_nlp_heurestic(x,y)

        x.upper_optimizer = x.upper_factory()
        upper_optimizer = x.upper_optimizer
        upper_vars = MOI.add_variables(upper_optimizer, x._variable_number)

        lower_variable_bounds = y.lower_variable_bounds
        upper_variable_bounds = y.upper_variable_bounds
        variable_info = x._input_problem._variable_info

        lvb = 0.0
        uvb = 0.0
        x0 = 0.0
        for i = 1:x._variable_number
            var = @inbounds variable_info[i]
            sv = MOI.SingleVariable(@inbounds upper_vars[i])
            if var.is_integer
            else
                lvb = @inbounds lower_variable_bounds[i]
                uvb = @inbounds upper_variable_bounds[i]
                if var.is_fixed
                    MOI.add_constraint(upper_optimizer, sv, ET(lvb))
                elseif var.has_lower_bound
                    if var.has_upper_bound
                        MOI.add_constraint(upper_optimizer, sv, LT(uvb))
                        MOI.add_constraint(upper_optimizer, sv, GT(lvb))
                    else
                        MOI.add_constraint(upper_optimizer, sv, GT(lvb))
                    end
                elseif var.has_upper_bound
                    MOI.add_constraint(upper_optimizer, sv, LT(uvb))
                end
                x0 = 0.5*(lvb + uvb)
                MOI.set(upper_optimizer, MOI.VariablePrimalStart(), @inbounds upper_vars[i], x0)
            end
        end

        # Add linear and quadratic constraints to model
        for (func, set) in x._input_problem._linear_leq_constraints
             MOI.add_constraint(upper_optimizer, func, set)
        end
        for (func, set) in x._input_problem._linear_geq_constraints
            MOI.add_constraint(upper_optimizer, func, set)
        end
        for (func, set) in x._input_problem._linear_eq_constraints
            MOI.add_constraint(upper_optimizer, func, set)
        end

        for (func, set) in x._input_problem._quadratic_leq_constraints
            MOI.add_constraint(upper_optimizer, func, set)
        end
        for (func, set) in x._input_problem._quadratic_geq_constraints
            MOI.add_constraint(upper_optimizer, func, set)
        end
        for (func, set) in x._input_problem._quadratic_eq_constraints
            MOI.add_constraint(upper_optimizer, func, set)
        end

        # Add nonlinear evaluation block
        MOI.set(upper_optimizer, MOI.NLPBlock(), x._input_problem._nlp_data)
        MOI.set(upper_optimizer, MOI.ObjectiveSense(), MOI.MIN_SENSE)

        if x._input_problem._objective_type === SINGLE_VARIABLE
            MOI.set(upper_optimizer, MOI.ObjectiveFunction{SV}(), x._input_problem._objective_sv)
        elseif x._input_problem._objective_type === SCALAR_AFFINE
            MOI.set(upper_optimizer, MOI.ObjectiveFunction{SAF}(), x._input_problem._objective_saf)
        elseif x._input_problem._objective_type === SCALAR_QUADRATIC
            MOI.set(upper_optimizer, MOI.ObjectiveFunction{SQF}(), x._input_problem._objective_sqf)
        end

        # Optimizes the object
        MOI.optimize!(upper_optimizer)

        # Process output info and save to CurrentUpperInfo object
        x._upper_termination_status = MOI.get(upper_optimizer, MOI.TerminationStatus())
        x._upper_result_status = MOI.get(upper_optimizer, MOI.PrimalStatus())

        if is_feasible_solution(x._upper_termination_status, x._upper_result_status)
            x._upper_feasibility = true
            value = MOI.get(upper_optimizer, MOI.ObjectiveValue())
            stored_adjusted_upper_bound!(x, value)
            x._best_upper_value = min(value, x._best_upper_value)
            x._upper_solution .= MOI.get(upper_optimizer, MOI.VariablePrimal(), upper_vars)
        else
            x._upper_feasibility = false
            x._upper_objective_value = Inf
        end
    else
        x._upper_feasibility = false
        x._upper_objective_value = Inf
    end
    return
end

"""
$(SIGNATURES)

Default upper bounding problem which simply calls `solve_local_nlp!` to solve
the nlp locally.
"""
function upper_problem!(t::ExtensionType, x::Optimizer)

    solve_local_nlp!(x)
    return
end


"""
$(SIGNATURES)

Default postprocess perfoms duality-based bound tightening on the `y`.
"""
function postprocess!(t::ExtensionType, x::Optimizer)

    if x._parameters.dbbt_depth > x._iteration_count
        variable_dbbt!(x._current_node, x._lower_lvd, x._lower_uvd,
                       x._lower_objective_value, x._global_upper_bound,
                       x._variable_number)
    end

    x._postprocess_feasibility = true

    return
end

"""
$(SIGNATURES)

Provides a hook for extensions to EAGO as opposed to standard global, local,
or linear solvers.
"""
function optimize_hook!(t::ExtensionType, x::Optimizer)
    return
end
