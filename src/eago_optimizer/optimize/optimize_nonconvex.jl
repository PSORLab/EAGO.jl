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
# TODO
#############################################################################

"""
$(TYPEDSIGNATURES)

Creates an initial node with initial box constraints and adds it to the stack.
"""
function create_initial_node!(m::Optimizer)

    branch_variable_count = m._branch_variable_count

    variable_info = m._working_problem._variable_info
    lower_bound = zeros(Float64, branch_variable_count)
    upper_bound = zeros(Float64, branch_variable_count)
    branch_count = 1

    for i = 1:m._working_problem._variable_count
        vi = @inbounds variable_info[i]
        if vi.branch_on === BRANCH
            @inbounds lower_bound[branch_count] = vi.lower_bound
            @inbounds upper_bound[branch_count] = vi.upper_bound
            branch_count += 1
        end
    end

    n = NodeBB(lower_bound, upper_bound, -Inf, Inf, 1, 1)
    push!(m._stack, n)
    m._node_count = 1
    m._maximum_node_id += 1

    return nothing
end

"""
$(TYPEDSIGNATURES)

Loads variables, linear constraints, and empty storage for first nlp and
quadratic cut.
"""
function load_relaxed_problem!(m::Optimizer)

    relaxed_optimizer = m.relaxed_optimizer

    # add variables and indices and constraints
    wp = m._working_problem
    variable_count = wp._variable_count
    branch_variable_count = 0

    for i = 1:variable_count

        relaxed_variable_indx = MOI.add_variable(relaxed_optimizer)
        relaxed_variable = SV(relaxed_variable_indx)
        push!(m._relaxed_variable_index, relaxed_variable_indx)

        vinfo = @inbounds wp._variable_info[i]

        is_branch_variable = @inbounds m._branch_variables[i]
        vinfo.branch_on = is_branch_variable ? BRANCH : NO_BRANCH
        is_branch_variable && (branch_variable_count += 1)

        if vinfo.is_integer

        elseif vinfo.is_fixed
            ci_sv_et = MOI.add_constraint(relaxed_optimizer, relaxed_variable, ET(vinfo.lower_bound))
            if is_branch_variable
                push!(m._relaxed_variable_eq, (ci_sv_et, branch_variable_count))
                wp._var_eq_count += 1
            end

        else
            if vinfo.has_lower_bound
                ci_sv_gt = MOI.add_constraint(relaxed_optimizer, relaxed_variable, GT(vinfo.lower_bound))
                if is_branch_variable
                    push!(m._relaxed_variable_gt, (ci_sv_gt, branch_variable_count))
                    wp._var_geq_count += 1
                end
            end

            if vinfo.has_upper_bound
                ci_sv_lt = MOI.add_constraint(relaxed_optimizer, relaxed_variable, LT(vinfo.upper_bound))
                if is_branch_variable
                    push!(m._relaxed_variable_lt, (ci_sv_lt, branch_variable_count))
                    wp._var_leq_count += 1
                end
            end
        end
    end

    m._branch_variable_count = branch_variable_count

    # add linear constraints
    add_linear_constraints!(m, relaxed_optimizer)

    MOI.set(relaxed_optimizer, MOI.ObjectiveSense(), MOI.MIN_SENSE)

    return nothing
end

function presolve_global!(t::ExtensionType, m::Optimizer)

    load_relaxed_problem!(m)
    create_initial_node!(m)

    branch_variable_count = m._branch_variable_count

    m._current_xref = fill(0.0, branch_variable_count)
    m._cut_solution = fill(0.0, branch_variable_count)
    m._lower_solution = fill(0.0, branch_variable_count)
    m._upper_solution = fill(0.0, branch_variable_count)
    m._lower_lvd = fill(0.0, branch_variable_count)
    m._lower_uvd = fill(0.0, branch_variable_count)
    m._continuous_solution = zeros(Float64, branch_variable_count)

    m._presolve_time = time() - m._parse_time

    return nothing
end

"""
$(SIGNATURES)

Selects node with the lowest lower bound in stack.
"""
function node_selection!(t::ExtensionType, m::Optimizer)
    m._node_count -= 1
    m._current_node = popmin!(m._stack)
    return nothing
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
    continue_flag = !isempty(m._stack)
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
  if (U < Inf) && (L > Inf)
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
        #for func in
        #end
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

    return nothing
end

"""
$(SIGNATURES)

Updates the relaxed constraint by setting the constraint set of `v == x*`` ,
`xL_i <= x_i`, and `x_i <= xU_i` for each such constraint added to the relaxed
optimizer.
"""
function update_relaxed_problem_box!(m::Optimizer, n::NodeBB)

    opt = m.relaxed_optimizer
    wp = m._working_problem

    lower_bound = n.lower_variable_bounds
    upper_bound = n.upper_variable_bounds

    relaxed_variable_eq = m._relaxed_variable_eq
    for i = 1:wp._var_eq_count
        constr_indx, node_indx = @inbounds relaxed_variable_eq[i]
        MOI.set(opt, MOI.ConstraintSet(), constr_indx, ET(@inbounds lower_bound[node_indx]))
    end

    relaxed_variable_lt = m._relaxed_variable_lt
    for i = 1:wp._var_leq_count
        constr_indx, node_indx = @inbounds relaxed_variable_lt[i]
        MOI.set(opt, MOI.ConstraintSet(), constr_indx, LT(@inbounds upper_bound[node_indx]))
    end

    relaxed_variable_gt = m._relaxed_variable_gt
    for i = 1:wp._var_geq_count
        constr_indx, node_indx = @inbounds relaxed_variable_gt[i]
        MOI.set(opt, MOI.ConstraintSet(), constr_indx, GT(@inbounds lower_bound[node_indx]))
    end

    return nothing
end

function interval_objective_bound(m::Optimizer, n::NodeBB)

    interval_objective_bound = bound_objective(m)

    if interval_objective_bound > m._lower_objective_value
        m._lower_objective_value = interval_objective_bound
        unsafe_check_fill!(m._lower_lvd, 0.0, m._relaxed_variable_number)
        unsafe_check_fill!(m._lower_uvd, 0.0, m._relaxed_variable_number)
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
        m._lower_objective_value = -Inf
    end
    m._lower_feasibility = feasible_flag

    return
end

"""
$(SIGNATURES)

Constructs and solves the relaxation using the default EAGO relaxation scheme
and optimizer on node `y`.
"""
function lower_problem!(t::ExtensionType, m::Optimizer)

    n = m._current_node

    if !m._obbt_performed_flag
        @__dot__ m._current_xref = 0.5*(n.lower_variable_bounds + n.upper_variable_bounds)
        unsafe_check_fill!(isnan, m._current_xref, 0.0, m._relaxed_variable_number)
        update_relaxed_problem_box!(m, n)
        relax_problem!(m, m._current_xref, 1)
    end
    relax_objective!(m, m._current_xref)

    # Optimizes the object
    relaxed_optimizer = m.relaxed_optimizer
    MOI.optimize!(relaxed_optimizer)

    m._lower_termination_status = MOI.get(relaxed_optimizer, MOI.TerminationStatus())
    m._lower_result_status = MOI.get(relaxed_optimizer, MOI.PrimalStatus())
    valid_flag, feasible_flag = is_globally_optimal(m._lower_termination_status, m._lower_result_status)

    if valid_flag && feasible_flag
        set_dual!(m)
        m._cut_add_flag = true
        m._lower_feasibility = true
        m._lower_objective_value = MOI.get(relaxed_optimizer, MOI.ObjectiveValue())
        for i = 1:m._relaxed_variable_number
            @inbounds m._lower_solution[i] = MOI.get(opt, MOI.VariablePrimal(), m._relaxed_variable_index[i])
        end

    elseif valid_flag
        m._cut_add_flag = false
        m._lower_feasibility  = false
        m._lower_objective_value = -Inf

    else
        fallback_interval_lower_bound!(m, n)
    end

    return nothing
end

"""
$(SIGNATURES)

Updates the internal storage in the optimizer after a valid feasible cut is added.
"""
function cut_update!(m::Optimizer)

    m._cut_feasibility = true

    relaxed_optimizer = x.relaxed_optimizer
    obj_val = MOI.get(relaxed_optimizer, MOI.ObjectiveValue())
    prior_obj_val = (x._cut_iterations == 2) ? x._lower_objective_value : x._cut_objective_value

    if prior_obj_val < obj_val

        x._cut_add_flag = true
        x._cut_objective_value = obj_val
        x._lower_objective_value = obj_val
        x._lower_termination_status = x._cut_termination_status
        x._lower_result_status = x._cut_result_status

        @inbounds x._cut_solution[:] = MOI.get(opt, MOI.VariablePrimal(), x._relaxed_variable_index)
        copyto!(x._lower_solution, x._cut_solution)
        set_dual!(x)
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
function cut_condition(t::ExtensionType, m::Optimizer)

    continue_cut_flag = m._cut_add_flag
    continue_cut_flag &= (m._cut_iterations < m.cut_max_iterations)
    n = m._current_node

    if continue_cut_flag
        xsol = (m._cut_iterations > 1) ? m._cut_solution : m._lower_solution
        xnew = (1.0 - m.cut_cvx)*mid(n) + m.cut_cvx*xsol
        if norm((m._current_xref - xnew)./diam(n), 1) > m.cut_tolerance
            m._current_xref = xnew
        else
            continue_cut_flag = false
        end
    end

    if !continue_cut_flag
        delete_nl_constraints!(m)
        delete_obj_cuts!(m)
    end

    # check to see if interval bound is preferable
    if m._lower_feasibility
        if m._objective_type === NONLINEAR
            objective_lo = eval_objective_lo(m._relaxed_evaluator)
        elseif m._objective_type === SINGLE_VARIABLE
                obj_indx = m._objective_sv.variable.value
                objective_lo = @inbounds y.lower_variable_bounds[obj_indx]
        elseif m._objective_type === SCALAR_AFFINE
                objective_lo = interval_bound(v._objective_saf, y, true)
        elseif m._objective_type === SCALAR_QUADRATIC
                objective_lo = interval_bound(m._objective_sqf, y, true)
        end
        if objective_lo > m._lower_objective_value
            m._lower_objective_value = objective_lo
            fill!(v._lower_lvd, 0.0)
            fill!(m._lower_uvd, 0.0)
        end
    end
    m._cut_iterations += 1

    return continue_cut_flag
end

"""
$(SIGNATURES)

Adds a cut for each constraint and the objective function to the subproblem.
"""
function add_cut!(t::ExtensionType, m::Optimizer)

    relax_problem!(m, m._current_xref, m._cut_iterations)
    relax_objective!(m, m._current_xref)

    # Optimizes the object
    relaxed_optimizer = m.relaxed_optimizer
    MOI.optimize!(opt)

    m._cut_termination_status = MOI.get(relaxed_optimizer, MOI.TerminationStatus())
    m._cut_result_status = MOI.get(relaxed_optimizer, MOI.PrimalStatus())
    valid_flag, feasible_flag = is_globally_optimal(m._cut_termination_status, m._cut_result_status)

    if valid_flag && feasible_flag
        cut_update!(m)
    elseif valid_flag
        m._cut_add_flag = false
        m._lower_feasibility  = false
        m._lower_objective_value = -Inf
    else
        m._cut_add_flag = false
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
function default_nlp_heurestic(m::Optimizer)
    bool = false
    ubd_limit = m._parameters.upper_bounding_depth
    depth = m._current_node.depth
    bool |= (depth <= ubd_limit)
    bool |= (rand() < 0.5^(depth - m._parameters.upper_bounding_depth))
    return bool
end

"""
$(SIGNATURES)

Default upper bounding problem which simply calls `solve_local_nlp!` to solve
the nlp locally.
"""
function upper_problem!(t::ExtensionType, m::Optimizer)

    if !default_nlp_heurestic(m)
        m._upper_feasibility = false
        m._upper_objective_value = Inf

    else
        solve_local_nlp!(m)

    end

    return nothing
end


"""
$(SIGNATURES)

Default postprocess perfoms duality-based bound tightening on the `y`.
"""
function postprocess!(t::ExtensionType, m::Optimizer)

    if m._parameters.dbbt_depth > m._iteration_count
        variable_dbbt!(m._current_node, m._lower_lvd, m._lower_uvd,
                       m._lower_objective_value, m._global_upper_bound,
                       m._variable_number)
    end

    return nothing
end

"""
$(SIGNATURES)

Provides a hook for extensions to EAGO as opposed to standard global, local,
or linear solvers.
"""
optimize_hook!(t::ExtensionType, m::Optimizer) = nothing

function store_candidate_solution!(m::Optimizer)
    if m._upper_feasibility && (m._upper_objective_value < m._global_upper_bound)
        m._feasible_solution_found = true
        m._first_solution_node = m._maximum_node_id
        m._solution_value = m._upper_objective_value
        m._global_upper_bound = m._upper_objective_value
        @__dot__ m._continuous_solution = m._upper_solution
    end
    return nothing
end

function set_global_lower_bound!(m::Optimizer)
    if !isempty(m._stack)
        min_node = minimum(m._stack)
        lower_bound = min_node.lower_bound
        if m._global_lower_bound < lower_bound
            m._global_lower_bound = lower_bound
        end
    end
    return nothing
end

# wraps subroutine call to isolate ExtensionType
parse_global!(m::Optimizer) = parse_global!(m.ext_type, m)
presolve_global!(m::Optimizer) = presolve_global!(m.ext_type, m)
termination_check(m::Optimizer) = termination_check(m.ext_type, m)
cut_condition(m::Optimizer) = cut_condition(m.ext_type, m)
convergence_check(m::Optimizer) = convergence_check(m.ext_type, m)
repeat_check(m::Optimizer) = repeat_check(m.ext_type, m)
node_selection!(m::Optimizer) = node_selection!(m.ext_type, m)
preprocess!(m::Optimizer) = preprocess!(m.ext_type, m)
lower_problem!(m::Optimizer) = lower_problem!(m.ext_type, m)
add_cut!(m::Optimizer) = add_cut!(m.ext_type, m)
upper_problem!(m::Optimizer) = upper_problem!(m.ext_type, m)
postprocess!(m::Optimizer) = postprocess!(m.ext_type, m)
single_storage!(m::Optimizer) = single_storage!(m.ext_type, m)
branch_node!(m::Optimizer) = branch_node!(m.ext_type, m)
fathom!(m::Optimizer) = fathom!(m.ext_type, m)

"""
$(TYPEDSIGNATURES)

Solves the branch and bound problem with the input EAGO optimizer object.
"""
function global_solve!(m::Optimizer)

    m._iteration_count = 1
    m._node_count = 1

    parse_global!(m)
    presolve_global!(m)

    logging_on = m._parameters.log_on
    verbosity = m._parameters.verbosity

    # terminates when max nodes or iteration is reach, or when node stack is empty
    while !termination_check(m)

        # Selects node, deletes it from stack, prints based on verbosity
        node_selection!(m)
        (verbosity >= 3) && print_node!(m)

        # Performs prepocessing and times
        logging_on && (start_time = time())
        preprocess!(m)
        if logging_on
            m._last_preprocess_time = time() - start_time
        end

        if m._preprocess_feasibility

            # solves & times lower bounding problem
            logging_on && (start_time = time())
            m._cut_iterations = 1
            lower_problem!(m)
            while cut_condition(m)
                add_cut!(m)
            end
            if logging_on
                m._last_lower_problem_time = time() - start_time
            end
            print_results!(m, true)
            print_results_post_cut!(m)

            # checks for infeasibility stores solution
            if m._lower_feasibility
                if !convergence_check(m)

                    logging_on && (start_time = time())
                    upper_problem!(m)
                    if logging_on
                        m._last_upper_problem_time = time() - start_time
                    end
                    print_results!(m, false)
                    store_candidate_solution!(m)
                    if m._input_problem._optimization_sense === MOI.FEASIBILITY_SENSE
                        if !m.feasible_local_continue || m.local_solve_only
                            break
                        end
                    end

                    # Performs and times post processing
                    logging_on && (start_time = time())
                    postprocess!(m)
                    if logging_on
                        m._last_postprocessing_time = time() - start_time
                    end

                    # Checks to see if the node
                    if m._postprocess_feasibility
                        if repeat_check(m)
                            single_storage!(m)
                        else
                            branch_node!(m)
                        end
                    end
                end
            end
            fathom!(m)
        else
            m._lower_objective_value = -Inf
            m._lower_feasibility = false
            m._upper_feasibility = false
        end
        set_global_lower_bound!(m)
        m._run_time = time() - m._start_time
        m._time_left = m.time_limit - m._run_time
        log_iteration!(m)
        print_iteration!(m)
        m._iteration_count += 1
    end

    m._objective_value = m._global_upper_bound

    # Prints the solution
    print_solution!(m)

    return nothing
end

optimize!(::Val{MINCVX}, m::Optimizer) = global_solve!(m)
