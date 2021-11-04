# Copyright (c) 2018: Matthew Wilhelm & Matthew Stuber.
# This code is licensed under MIT license (see LICENSE.md for full details)
#############################################################################
# EAGO
# A development environment for robust and global optimization
# See https://github.com/PSORLab/EAGO.jl
#############################################################################
# src/eago_optimizer/domain_reduction.jl
# Contains subroutines used for domain reduction.
#############################################################################

"""
$(FUNCTIONNAME)

Tightens the bounds of the `_current_node` using the current global upper bound
and the duality information obtained from the relaxation.
"""
function variable_dbbt!(n::NodeBB, mult_lo::Vector{Float64}, mult_hi::Vector{Float64},
                        LBD::Float64, UBD::Float64, nx::Int64)

    delta = UBD - LBD
    lvbs = n.lower_variable_bounds
    uvbs = n.upper_variable_bounds
    if LBD <= UBD
        for i = 1:nx
            ml = @inbounds mult_lo[i]
            if ml > 0.0
                cut = @inbounds lvbs[i] + delta/ml
                if cut < @inbounds uvbs[i]
                    @inbounds uvbs[i] = cut
                end
            else
                mh = @inbounds mult_hi[i]
                if mh > 0.0
                    cut = @inbounds uvbs[i] - delta/mh
                    if cut > @inbounds lvbs[i]
                        @inbounds lvbs[i] = cut
                    end
                end
            end
         end
    end

    return nothing
end

function set_preprocess_status(m::GlobalOptimizer{R,S,Q}, d) where {R,S,Q<:ExtensionType}
    tstatus = MOI.get(d, MOI.TerminationStatus())
    pstatus = MOI.get(d, MOI.PrimalStatus())
    dstatus =  MOI.get(d, MOI.DualStatus())
    m._preprocess_termination_status = tstatus
    m._preprocess_primal_status = pstatus
    m._preprocess_dual_status = dstatus
    return relaxed_problem_status(tstatus, pstatus, dstatus)
end

"""
$(FUNCTIONNAME)

Excludes OBBT on variable indices that are tight for the solution of the relaxation.
"""
function trivial_filtering!(m::GlobalOptimizer{R,S,Q}, n::NodeBB) where {R,S,Q<:ExtensionType}
    
    d = _relaxed_optimizer(m)
    obbt_tolerance = _obbt_tolerance(m)
    if set_preprocess_status(m,d) == RRS_OPTIMAL
        for j = 1:m._obbt_variable_count
            vi = VI(_bvi(m, m._relaxed_variable_index[j].value))
            if m._obbt_working_lower_index[j]
                z = MOI.get(d, MOI.VariablePrimal(),vi) - n.lower_variable_bounds[j]
                if abs(z) <= obbt_tolerance
                    m._obbt_working_lower_index[j] = false
                end
            end
            if m._obbt_working_upper_index[j]
                z = n.upper_variable_bounds[j] - MOI.get(d, MOI.VariablePrimal(), vi)
                if abs(z) <= obbt_tolerance
                    m._obbt_working_upper_index[j] = false
                end
            end
        end
    end
    return
end

and_not(x,y) = x & ~y
"""
$(FUNCTIONNAME)

Utility function used to set vector of booleans z to x & ~y. Avoids the
generation of conversion of the BitArray created by broadcasting logical operators.
"""
bool_indx_diff!(z::Vector{Bool},x::Vector{Bool}, y::Vector{Bool}) = map!(and_not, z, x, y)

"""
$(FUNCTIONNAME)

Excludes OBBT on variable indices after a search in a filtering direction.
"""
function aggressive_filtering!(m::GlobalOptimizer{R,S,Q}, n::NodeBB) where {R,S,Q<:ExtensionType}
    if _obbt_aggressive_on(m)
        # Initial filtering vector (negative one direction per remark in Gleixner2017)
        d = _relaxed_optimizer(m)
        variable_number = _variable_num(FullVar(), m)
        v = -ones(variable_number)

        # Copy prior index set (ignores linear and binary terms)
        obbt_variable_count = m._obbt_variable_count
        copyto!(m._old_low_index, m._obbt_working_lower_index)
        copyto!(m._old_upp_index, m._obbt_working_upper_index)
        copyto!(m._new_low_index, m._obbt_working_lower_index)
        copyto!(m._new_upp_index, m._obbt_working_upper_index)

        # Exclude unbounded directions
        @__dot__ m._new_low_index = m._new_low_index & !isinf(n.lower_variable_bounds)

        # Begin the main algorithm
        for k = 1:_obbt_aggressive_max_iteration(m)

            # Set index differences and vector for filtering direction
            @__dot__ m._lower_indx_diff = m._old_low_index & ~m._new_low_index
            @__dot__ m._upper_indx_diff = m._old_upp_index & ~m._new_upp_index

            for i = 1:obbt_variable_count
                vi = v[i]
                if m._lower_indx_diff[i] && vi < 0.0
                    v[i] = 0.0
                end
                if m._upper_indx_diff[i] && vi > 0.0
                    v[i] = 0.0
                end
            end

            # Termination Condition
            ((~any(m._new_low_index) & ~any(m._new_upp_index)) || (iszero(v))) && break
            if k >= 2
                if (count(m._lower_indx_diff) + count(m._upper_indx_diff)) < m._parameters.obbt_aggressive_min_dimension
                    break
                end
            end

            # Set objective in OBBT problem to filtering vector
            MOI.set(d, MOI.ObjectiveSense(), MOI.MAX_SENSE)
            MOI.set(d, MOI.ObjectiveFunction{SAF}(), SAF(SAT.(v, m._relaxed_variable_index), 0.0))

            # Optimizes the problem and if successful filter additional bounds
            MOI.optimize!(d)

            if set_preprocess_status(m,d) == RRS_OPTIMAL
                variable_primal = MOI.get(d, MOI.VariablePrimal(), m._relaxed_variable_index)
                copyto!(m._new_low_index, m._old_low_index)
                copyto!(m._new_upp_index, m._old_upp_index)
                for i = 1:obbt_variable_count
                    vp_value = variable_primal[i]
                    if m._old_low_index[i] && vp_value == n.lower_variable_bounds[i]
                        m._new_low_index[i] = false
                    end
                    if m._old_upp_index[i] && vp_value == n.upper_variable_bounds[i]
                        m._new_upp_index[i] = false
                    end
                end
            else
                return false
            end
        end
        copyto!(m._obbt_working_lower_index, m._new_low_index)
        copyto!(m._obbt_working_upper_index, m._new_upp_index)
    end
    return true
end

"""
$(FUNCTIONNAME)
"""
function set_reference_point!(m::GlobalOptimizer)

    wp = m._working_problem
    evaluator = wp._relaxed_evaluator
    evaluator_x = evaluator.variable_values.x
    current_xref = m._lower_solution

    new_reference_point = false
    for i = 1:length(current_xref)
        node_x = current_xref[i]
        solution_x = evaluator_x[i]
        if node_x !== solution_x
            evaluator_x[i] = node_x
            new_reference_point = true
        end
    end
    #=
    for node_i = 1:m._branch_variable_count
        solution_i = m._branch_to_sol_map[node_i]

        node_x = current_xref[node_i]
        solution_x = evaluator_x[solution_i]

        if node_x !== solution_x
            evaluator_x[solution_i] = node_x
            new_reference_point = true
        end
    end
    =#

    if new_reference_point
        foreach(c -> _set_has_value!(c, false), wp._nonlinear_constr)
        _set_has_value!(wp._objective, false)
    end
    fill!(evaluator.subexpressions_eval, false)

    return nothing
end

"""
$(FUNCTIONNAME)

Performs OBBT with filtering and greedy ordering as detailed in:
Gleixner, A.M., Berthold, T., Müller, B. et al. J Glob Optim (2017) 67: 731.
https://doi.org/10.1007/s10898-016-0450-4
"""
function obbt!(m::GlobalOptimizer{R,S,Q}) where {R,S,Q<:ExtensionType}

    feasibility = true
    n = m._current_node
    d = _relaxed_optimizer(m)

    if relax_problem!(m)
        MOI.set(d, MOI.ObjectiveSense(), MOI.MIN_SENSE)
        MOI.optimize!(d)

        # Sets indices to attempt OBBT on
        obbt_variable_count = m._obbt_variable_count
        fill!(m._obbt_working_lower_index, true)
        fill!(m._obbt_working_upper_index, true)

        # Filters out any indicies with active bounds on variables
        # determined by solving the feasibility problem
        trivial_filtering!(m, n)
        feasibility = aggressive_filtering!(m, n)
        if set_preprocess_status(m,d) == RRS_OPTIMAL
            xLP = MOI.get(d, MOI.VariablePrimal(), m._relaxed_variable_index)
        else
            return false
        end

        # continue tightening bounds by optimization until all indices have been checked
        # or the node is empty and the problem is thus proven infeasible
        while (any(m._obbt_working_lower_index) || any(m._obbt_working_upper_index)) && !isempty(n)

            # Get lower value
            lower_indx = -1;     upper_indx = -1
            lower_value = Inf;   upper_value = Inf

            # min of xLP - yL on active
            if any(m._obbt_working_lower_index)
                for i = 1:obbt_variable_count
                    if m._obbt_working_lower_index[i]
                        temp_value = _lower_solution(BranchVar(),m,i) - _lower_bound(BranchVar(),m,i)
                        if temp_value <= lower_value   # Need less than or equal to handle unbounded cases
                            lower_value = temp_value
                            lower_indx = i
                        end
                    end
                end
            end

            # min of yU - xLP on active
            if any(m._obbt_working_upper_index)
                for i = 1:obbt_variable_count
                    if m._obbt_working_upper_index[i]
                        temp_value = _upper_bound(BranchVar(),m,i) - _lower_solution(BranchVar(),m,i)
                        if temp_value <= upper_value
                            upper_value = temp_value
                            upper_indx = i
                        end
                    end
                end
            end

            # default to upper bound if no lower bound is found, use maximum distance otherwise
            if lower_value <= upper_value && lower_indx > 0
                m._obbt_working_lower_index[lower_indx] = false
                MOI.set(d, MOI.ObjectiveSense(), MOI.MIN_SENSE)
                MOI.set(d, MOI.ObjectiveFunction{SV}(), SV(m._relaxed_variable_index[_bvi(m, lower_indx)]))
                MOI.optimize!(d)

                status = set_preprocess_status(m,d)
                if status == RRS_OPTIMAL
                    xLP .= MOI.get(d, MOI.VariablePrimal(), m._relaxed_variable_index)
                    updated_value = xLP[_bvi(m, lower_indx)]
                    previous_value = n.lower_variable_bounds[lower_indx]

                    # if bound is improved update node and corresponding constraint update
                    # the node bounds and the single variable bound in the relaxation
                    # we assume branching does not occur on fixed variables and interval
                    # constraints are internally bridged by EAGO. So the only L <= x
                    # constraint in the model is a GreaterThan.
                    if updated_value > previous_value && (updated_value - previous_value) > 1E-6
                        sv_geq_ci = m._node_to_sv_geq_ci[lower_indx]
                        MOI.set(d, MOI.ConstraintSet(), sv_geq_ci, GT(updated_value))
                        n.lower_variable_bounds[lower_indx] = updated_value
                    end

                    if isempty(n)
                        feasibility = false
                        break
                    end

                elseif status == RRS_INFEASIBLE
                    feasibility = false
                    break
                else
                    break
                end

            elseif upper_indx > 0
                
                m._obbt_working_upper_index[upper_indx] = false
                MOI.set(d, MOI.ObjectiveSense(), MOI.MAX_SENSE)
                MOI.set(d, MOI.ObjectiveFunction{SV}(), SV(m._relaxed_variable_index[_bvi(m, upper_indx)]))
                MOI.optimize!(d)

                status = set_preprocess_status(m,d)
                if status == RRS_OPTIMAL
                    xLP .= MOI.get(d, MOI.VariablePrimal(), m._relaxed_variable_index)
                    updated_value = xLP[_bvi(m, upper_indx)]
                    previous_value = n.upper_variable_bounds[upper_indx]

                    # if bound is improved update node and corresponding constraint update
                    # the node bounds and the single variable bound in the relaxation
                    # we assume branching does not occur on fixed variables and interval
                    # constraints are internally bridged by EAGO. So the only U => x
                    # constraint in the model is a LessThan.
                    if updated_value < previous_value && (previous_value - updated_value) > 1E-6
                        sv_leq_ci = m._node_to_sv_leq_ci[upper_indx]
                        MOI.set(d, MOI.ConstraintSet(), sv_leq_ci, LT(updated_value))
                        n.upper_variable_bounds[upper_indx] = updated_value
                    end

                    if isempty(n)
                        feasibility = false
                        break
                    end

                elseif status == RRS_INFEASIBLE
                    feasibility = false
                    break
                else
                    break
                end
            else
                break
            end
            trivial_filtering!(m, n)
        end
    end
    return feasibility
end

"""
$(FUNCTIONNAME)
"""
function load_fbbt_buffer!(m::GlobalOptimizer)
    for i = 1:m._working_problem._variable_count
        m._lower_fbbt_buffer[i] = _lower_bound(FullVar(), m, i)
        m._upper_fbbt_buffer[i] = _upper_bound(FullVar(), m, i)
    end
    return
end

"""
$(FUNCTIONNAME)
"""
function unpack_fbbt_buffer!(m::GlobalOptimizer)

    n = m._current_node
    sol_to_branch = m._sol_to_branch_map
    lower_variable_bounds = n.lower_variable_bounds
    upper_variable_bounds = n.upper_variable_bounds

    for i = 1:m._working_problem._variable_count
        if m._branch_variables[i]
            indx = sol_to_branch[i]
            if m._lower_fbbt_buffer[i] > lower_variable_bounds[indx]
                lower_variable_bounds[indx] = m._lower_fbbt_buffer[i]
            end
            if upper_variable_bounds[indx] > m._upper_fbbt_buffer[i]
                upper_variable_bounds[indx] = m._upper_fbbt_buffer[i]
            end
        end
    end

    return nothing
end

"""
$(FUNCTIONNAME)

Performs feasibility-based bound tightening on a back-end constraint and returns `true` if it is feasible or
`false` if it is infeasible.
"""
function fbbt! end

function fbbt!(m::GlobalOptimizer, f::AffineFunctionIneq)

    # compute full sum
    lower_bounds = m._lower_fbbt_buffer
    upper_bounds = m._upper_fbbt_buffer

    terms = f.terms
    temp_sum = -f.constant

    for k = 1:f.len

        aik, indx_k = @inbounds terms[k]
        if aik !== 0.0
            aik_xL = aik*(@inbounds lower_bounds[indx_k])
            aik_xU = aik*(@inbounds upper_bounds[indx_k])
            temp_sum -= min(aik_xL, aik_xU)
        end
    end

    # subtract extra term, check to see if implied bound is better, if so update the node and
    # the working sum if the node is now empty then break

    for k = 1:f.len

        aik, indx_k = @inbounds terms[k]
        if aik !== 0.0

            xL = @inbounds lower_bounds[indx_k]
            xU = @inbounds upper_bounds[indx_k]

            aik_xL = aik*xL
            aik_xU = aik*xU

            temp_sum += min(aik_xL, aik_xU)
            xh = temp_sum/aik

            if aik > 0.0
                (xh < xL) && return false
                if xh > xL
                    @inbounds upper_bounds[indx_k] = xh
                end

            elseif aik < 0.0
                (xh > xU) && return false
                if xh < xU
                    @inbounds lower_bounds[indx_k] = xh
                end

            else
                temp_sum -= min(aik_xL, aik_xU)
                continue

            end

            aik_xL = aik*(@inbounds lower_bounds[indx_k])
            aik_xU = aik*(@inbounds upper_bounds[indx_k])
            temp_sum -= min(aik_xL, aik_xU)

        end
    end

    return true
end

function fbbt!(m::GlobalOptimizer, f::AffineFunctionEq)

    # compute full sum
    lower_bounds = m._lower_fbbt_buffer
    upper_bounds = m._upper_fbbt_buffer

    terms = f.terms
    temp_sum_leq = -f.constant
    temp_sum_geq = -f.constant

    for k = 1:f.len
        aik, indx_k = @inbounds terms[k]

        if aik !== 0.0
            aik_xL = aik*(@inbounds lower_bounds[indx_k])
            aik_xU = aik*(@inbounds upper_bounds[indx_k])
            temp_sum_leq -= min(aik_xL, aik_xU)
            temp_sum_geq -= max(aik_xL, aik_xU)

        end
    end

    # subtract extra term, check to see if implied bound is better, if so update the node and
    # the working sum if the node is now empty then break
    for k = 1:f.len

        aik, indx_k = @inbounds terms[k]
        if aik !== 0.0

            xL = @inbounds lower_bounds[indx_k]
            xU = @inbounds upper_bounds[indx_k]

            aik_xL = aik*xL
            aik_xU = aik*xU

            temp_sum_leq += min(aik_xL, aik_xU)
            temp_sum_geq += max(aik_xL, aik_xU)

            xh_leq = temp_sum_leq/aik
            xh_geq = temp_sum_geq/aik

            if aik > 0.0
                (xh_leq < xL) && return false
                if xh_leq > xL
                    @inbounds upper_bounds[indx_k] = xh_leq
                end
                (xh_geq > xU) && return false
                if (xh_geq < xU)
                    @inbounds lower_bounds[indx_k] = xh_geq
                end

            elseif aik < 0.0
                (xh_leq > xU) && return false
                if xh_leq < xU
                    @inbounds lower_bounds[indx_k] = xh_leq
                end
                (xh_geq < xL) && return false
                if (xh_geq > xL)
                    @inbounds upper_bounds[indx_k] = xh_geq
                end

            else
                temp_sum_leq -= min(aik_xL, aik_xU)
                temp_sum_geq -= max(aik_xL, aik_xU)
                continue

            end
            aik_xL = aik*(@inbounds lower_bounds[indx_k])
            aik_xU = aik*(@inbounds upper_bounds[indx_k])

            temp_sum_leq -= min(aik_xL, aik_xU)
            temp_sum_geq -= max(aik_xL, aik_xU)
        end
    end

    return true
end

cp_condition(m::GlobalOptimizer) = false

_propagate_constraint!(d, f) = true
function _propagate_constraint!(d, f::BufferedNonlinearFunction)
    forward_pass!(d, f)
    is_feasible = rprop!(Relax(), d, f)
    d.interval_intersect = true
    is_feasible && forward_pass!(d, f)
end
"""
Performs bound tightening based on forward/reverse interval and/or McCormick passes. This routine
resets the current node with new interval bounds.
"""
function set_constraint_propagation_fbbt!(m::GlobalOptimizer{R,S,Q}) where {R,S,Q<:ExtensionType}
    feasible_flag = true

    wp = m._working_problem
    if m._nonlinear_evaluator_created
        evaluator = wp._relaxed_evaluator
        set_node!(wp._relaxed_evaluator, m._current_node)
        set_reference_point!(m)

        wp._relaxed_evaluator.is_first_eval = m._new_eval_constraint
        for constr in wp._nonlinear_constr
            if feasible_flag
                forward_pass!(evaluator, constr)
                feasible_flag &= rprop!(Relax(), evaluator, constr)
                evaluator.interval_intersect = true
            end
        end

        wp._relaxed_evaluator.is_first_eval = m._new_eval_constraint
        for constr in wp._nonlinear_constr
            feasible_flag && forward_pass!(evaluator, constr)
        end

        evaluator.is_post = m._parameters.subgrad_tighten
        wp._relaxed_evaluator.is_first_eval = m._new_eval_objective
        feasible_flag && _propagate_constraint!(evaluator, wp._objective)

        m._new_eval_constraint = false
        m._new_eval_objective = false
        _get_x!(BranchVar, m._current_xref, evaluator)
        m._current_node = retrieve_node(evaluator)
    end

    return feasible_flag
end
