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
    return
end

"""
$(FUNCTIONNAME)

Excludes OBBT on variable indices that are tight for the solution of the relaxation.
"""
function trivial_filtering!(m::Optimizer, n::NodeBB)

    obbt_tolerance = m._parameters.obbt_tolerance
    m._preprocess_termination_status = MOI.get(m.relaxed_optimizer, MOI.TerminationStatus())
    m._preprocess_result_status = MOI.get(m.relaxed_optimizer, MOI.PrimalStatus())
    valid_flag, feasible_flag = is_globally_optimal(m._preprocess_termination_status,
                                                    m._preprocess_result_status)

    if valid_flag && feasible_flag
        for j = 1:length(m._obbt_working_lower_index)
            if @inbounds m._obbt_working_lower_index[j]
                vi = @inbounds m._relaxed_variable_index[j]
                diff = MOI.get(m.relaxed_optimizer, MOI.VariablePrimal(), vi)
                diff -= @inbounds n.lower_variable_bounds[j]
                if abs(diff) <= obbt_tolerance
                    @inbounds m._obbt_working_lower_index[j] = false
                end
            end
        end
        for j = 1:length(m._obbt_working_upper_index)
            if @inbounds m._obbt_working_upper_index[j]
                vi = @inbounds m._relaxed_variable_index[j]
                diff = -MOI.get(m.relaxed_optimizer, MOI.VariablePrimal(), vi)
                diff += @inbounds n.upper_variable_bounds[j]
                if abs(diff) <= obbt_tolerance
                    @inbounds m._obbt_working_upper_index[j] = false
                end
            end
        end
    end

    return
end

"""
$(FUNCTIONNAME)

Utility function used to set vector of booleans z to x & ~y. Avoids the
generation of conversion of the BitArray created by broadcasting logical operators.
"""
function bool_indx_diff(z::Vector{Bool},x::Vector{Bool}, y::Vector{Bool})
    for i = 1:length(z)
        @inbounds z[i] = (x[i] & ~y[i])
    end
    return
end

"""
$(FUNCTIONNAME)

Excludes OBBT on variable indices after a search in a filtering direction.
"""
function aggressive_filtering!(m::Optimizer, n::NodeBB)

    # Initial filtering vector (negative one direction per remark in Gleixner2017)
    variable_number = m._working_problem._variable_count
    v = -ones(variable_number)

    # Copy prior index set (ignores linear and binary terms)
    obbt_variable_count = m._obbt_variable_count
    copyto!(m._old_low_index, m._obbt_working_lower_index)
    copyto!(m._old_upp_index, m._obbt_working_upper_index)
    copyto!(m._new_low_index, m._obbt_working_lower_index)
    copyto!(m._new_upp_index, m._obbt_working_upper_index)

    # Exclude unbounded directions
    for i = 1:obbt_variable_count
        if @inbounds m._new_low_index[i] && @inbounds n.lower_variable_bounds[i] === -Inf
            @inbounds m._new_low_index[i] = false
        end
        if @inbounds m._new_low_index[i] && @inbounds n.upper_variable_bounds[i] === Inf
            @inbounds m._new_low_index[i] = false
        end
    end

    # Begin the main algorithm
    for k = 1:m._parameters.obbt_aggressive_max_iteration

        # Set index differences and vector for filtering direction
        bool_indx_diff(m._lower_indx_diff, m._old_low_index, m._new_low_index)
        bool_indx_diff(m._upper_indx_diff, m._old_upp_index, m._new_upp_index)

        for i = 1:obbt_variable_count
            vi = @inbounds v[i]
            if @inbounds m._lower_indx_diff[i] && vi < 0.0
                @inbounds v[i] = 0.0
            end
            if @inbounds m._upper_indx_diff[i] && vi > 0.0
                @inbounds v[i] = 0.0
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
        MOI.set(m.relaxed_optimizer, MOI.ObjectiveSense(), MOI.MAX_SENSE)
        saf = SAF(SAT.(v, m._relaxed_variable_index), 0.0)
        MOI.set(m.relaxed_optimizer, MOI.ObjectiveFunction{SAF}(), saf)

        # Optimizes the problem and if successful filter additional bounds
        MOI.optimize!(m.relaxed_optimizer)

        m._preprocess_termination_status = MOI.get(m.relaxed_optimizer, MOI.TerminationStatus())
        m._preprocess_result_status = MOI.get(m.relaxed_optimizer, MOI.PrimalStatus())
        valid_flag, feasible_flag = is_globally_optimal(m._preprocess_termination_status,
                                                        m._preprocess_result_status)

        if valid_flag && feasible_flag
            variable_primal = MOI.get(m.relaxed_optimizer, MOI.VariablePrimal(), m._relaxed_variable_index)
            copyto!(m._new_low_index, m._old_low_index)
            copyto!(m._new_upp_index, m._old_upp_index)
            for i = 1:obbt_variable_count
                vp_value =  @inbounds variable_primal[i]
                if @inbounds m._old_low_index[i] && vp_value == @inbounds n.lower_variable_bounds[i]
                    @inbounds m._new_low_index[i] = false
                end
                if @inbounds m._old_upp_index[i] && vp_value == @inbounds n.upper_variable_bounds[i]
                    @inbounds m._new_upp_index[i] = false
                end
            end
        else
            return false
        end
    end
    copyto!(m._obbt_working_lower_index, m._new_low_index)
    copyto!(m._obbt_working_upper_index, m._new_upp_index)
    return true
end

"""
$(FUNCTIONNAME)

Performs OBBT with filtering and greedy ordering as detailed in:
Gleixner, A.M., Berthold, T., Müller, B. et al. J Glob Optim (2017) 67: 731.
https://doi.org/10.1007/s10898-016-0450-4
"""
function obbt!(m::Optimizer)

    feasibility = true

    n = m._current_node
    branch_to_sol_map = m._branch_to_sol_map
    relaxed_optimizer = m.relaxed_optimizer

    # solve initial problem to feasibility
    set_first_relax_point!(m)
    update_relaxed_problem_box!(m)
    relax_constraints!(m, 1)
    relax_objective!(m, 1)

    MOI.set(relaxed_optimizer, MOI.ObjectiveSense(), MOI.FEASIBILITY_SENSE)
    MOI.optimize!(relaxed_optimizer)

    # Sets indices to attempt OBBT on (full set...)
    obbt_variable_count = m._obbt_variable_count
    obbt_variables = m._obbt_variables
    fill!(m._obbt_working_lower_index, true)
    fill!(m._obbt_working_upper_index, true)

    # Prefiltering steps && and sets initial LP values
    trivial_filtering!(m, n)
    if m._parameters.obbt_aggressive_on
        feasibility = aggressive_filtering!(m, n)
    end

    m._preprocess_termination_status = MOI.get(relaxed_optimizer, MOI.TerminationStatus())
    m._preprocess_result_status = MOI.get(relaxed_optimizer, MOI.PrimalStatus())
    valid_flag, feasible_flag = is_globally_optimal(m._preprocess_termination_status,
                                                    m._preprocess_result_status)

    if valid_flag && feasible_flag
        xLP = MOI.get(relaxed_optimizer, MOI.VariablePrimal(), m._relaxed_variable_index)
    else
        return false
    end

    while (any(m._obbt_working_lower_index) || any(m._obbt_working_upper_index)) & ~isempty(n)

        # Get lower value
        lower_indx = -1;     upper_indx = -1
        lower_value = Inf;   upper_value = Inf

        # min of xLP - yL on active
        if any(m._obbt_working_lower_index)
            for i = 1:obbt_variable_count
                if @inbounds m._obbt_working_lower_index[i]
                    sol_indx = branch_to_sol_map[i]
                    temp_value = @inbounds xLP[sol_indx] - n.lower_variable_bounds[i]
                    # Need less than or equal to handle unbounded cases
                    if temp_value <= lower_value
                        lower_value = temp_value
                        lower_indx = i
                    end
                end
            end
        end

        # min of yU - xLP on active
        if any(m._obbt_working_upper_index)
            for i = 1:obbt_variable_count
                if @inbounds m._obbt_working_upper_index[i]
                    sol_indx = branch_to_sol_map[i]
                    temp_value = @inbounds n.upper_variable_bounds[i] - xLP[sol_indx]
                    if temp_value <= upper_value
                        upper_value = temp_value
                        upper_indx = i
                    end
                end
            end
        end

        # default to upper bound if no lower bound is found, use maximum distance otherwise
        if lower_value <= upper_value && lower_indx > 0

            @inbounds m._obbt_working_lower_index[lower_indx] = false
            var = SV(m._relaxed_variable_index[lower_indx])

            MOI.set(relaxed_optimizer, MOI.ObjectiveSense(), MOI.MIN_SENSE)
            MOI.set(relaxed_optimizer, MOI.ObjectiveFunction{SV}(), var)

            MOI.optimize!(m.relaxed_optimizer)
            m._preprocess_termination_status = MOI.get(relaxed_optimizer, MOI.TerminationStatus())
            m._preprocess_result_status = MOI.get(relaxed_optimizer, MOI.PrimalStatus())
            valid_flag, feasible_flag = is_globally_optimal(m._preprocess_termination_status,
                                                            m._preprocess_result_status)

            if valid_flag && feasible_flag
                xLP .= MOI.get(relaxed_optimizer, MOI.VariablePrimal(), m._relaxed_variable_index)
                node_index = branch_to_sol_map[lower_indx]
                if m._working_problem._variable_info[node_index].is_integer
                    @inbounds n.lower_variable_bounds[lower_indx] = ceil(xLP[node_index])
                else
                    @inbounds n.lower_variable_bounds[lower_indx] = xLP[node_index]
                end
                if isempty(n)
                    feasibility = false
                    break
                end
            elseif valid_flag
                feasibility = false
            else
                break
            end

        elseif upper_indx > 0

            m._obbt_working_upper_index[upper_indx] = false
            var = SV(m._relaxed_variable_index[upper_indx])
            MOI.set(relaxed_optimizer, MOI.ObjectiveSense(), MOI.MAX_SENSE)
            MOI.set(relaxed_optimizer, MOI.ObjectiveFunction{SV}(), var)
            MOI.optimize!(relaxed_optimizer)
            m._preprocess_termination_status = MOI.get(m.relaxed_optimizer, MOI.TerminationStatus())
            m._preprocess_result_status = MOI.get(m.relaxed_optimizer, MOI.PrimalStatus())
            valid_flag, feasible_flag = is_globally_optimal(m._preprocess_termination_status,
                                                            m._preprocess_result_status)

            if valid_flag && feasible_flag
                xLP .= MOI.get(relaxed_optimizer, MOI.VariablePrimal(), m._relaxed_variable_index)
                node_index = branch_to_sol_map[upper_indx]
                if m._working_problem._variable_info[node_index].is_integer
                    @inbounds n.upper_variable_bounds[upper_indx] = ceil(xLP[node_index])
                else
                    @inbounds n.upper_variable_bounds[upper_indx] = xLP[node_index]
                end
                if isempty(n)
                    feasibility = false
                    break
                end
            elseif valid_flag
                feasibility = false
            else
                break
            end

        else
            break
        end
        trivial_filtering!(m, n)
    end

    return feasibility
end


function load_fbbt_buffer!(m::Optimizer)
end

function unpack_fbbt_buffer!(m::Optimizer)
end

"""

Performs feasibility-based bound tightening on a back-end constraint and returns `true` if it is feasible or
`false` if it is infeasible.
"""
function fbbt! end

function fbbt!(m::Optimizer, f::AffineFunctionIneq)

    println("ran me affine ineq")
    # compute full sum
    n = m._current_node
    lower_bounds = n._lower_fbbt_buffer
    upper_bounds = n._upper_fbbt_buffer

    terms = f.terms
    temp_sum = f.constant

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

            aik_xL = aik*(@inbounds lower_bounds[indx_k])
            aik_xU = aik*(@inbounds upper_bounds[indx_k])
            temp_sum += min(aik_xL, aik_xU)
            xh = temp_sum/aik

            if aik > 0.0
                (xh > xU) && return false
                if xh > xL
                    @inbounds lower_bounds[indx_k] = xh
                end

            elseif aik < 0.0
                (xh < xL) && return false
                if xh < xU
                    @inbounds upper_bounds[indx_k] = xh
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

function fbbt!(m::Optimizer, f::AffineFunctionEq)

    # compute full sum
    n = m._current_node
    lower_bounds = n._lower_fbbt_buffer
    upper_bounds = n._upper_fbbt_buffer

    terms = f.terms
    temp_sum_leq = f.constant
    temp_sum_geq = -f.constant

    for k = 1:f.len

        aik, indx_k = @inbounds terms[k]
        if aik !== 0.0
            aik_xL = aik*(@inbounds lower_bounds[indx_k])
            aik_xU = aik*(@inbounds upper_bounds[indx_k])
            temp_sum_leq -= min(aik_xL, aik_xU)
            temp_sum_geq += max(aik_xL, aik_xU)
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
            xh_leq = temp_sum_leq/aik

            if aik > 0.0
                (xh_leq > xU) && return false
                if xh_leq > xL
                    @inbounds lower_bounds[indx_k] = xh_leq
                    aik_xL = aik*xh_leq
                end
                temp_sum_geq -= max(aik_xL, aik_xU)
                xh_geq = -temp_sum_geq/aik
                (xh_geq < xL) && return false
                (xh_geq > xU) && (@inbounds lower_bounds[indx_k] = xh_leq)

            elseif aik < 0.0
                (xh_leq < xL) && return false
                if xh_leq < xU
                    @inbounds upper_bounds[indx_k] = xh_leq
                    aik_xU = aik*xh_leq
                end
                temp_sum_geq -= max(aik_xL, aik_xU)
                xh_geq = -temp_sum_geq/aik
                (xh_geq > xU) && return false
                (xh_geq > xL) && (@inbounds lower_bounds[indx_k] = xh_geq)

            else
                temp_sum_leq -= min(aik_xL, aik_xU)
                temp_sum_geq += max(aik_xL, aik_xU)
                continue

            end
            aik_xL = aik*(@inbounds lower_bounds[indx_k])
            aik_xU = aik*(@inbounds upper_bounds[indx_k])
            temp_sum_leq -= min(aik_xL, aik_xU)
            temp_sum_geq += max(aik_xL, aik_xU)
        end
    end

    return true
end

cp_condition(m::Optimizer) = false

function set_constraint_propagation_fbbt(m::Optimizer)
end
