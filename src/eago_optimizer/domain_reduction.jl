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
function trivial_filtering!(x::Optimizer, y::NodeBB)

    x._preprocess_termination_status = MOI.get(x.relaxed_optimizer, MOI.TerminationStatus())
    x._preprocess_result_status = MOI.get(x.relaxed_optimizer, MOI.PrimalStatus())
    valid_flag, feasible_flag = is_globally_optimal(x._preprocess_termination_status,
                                                    x._preprocess_result_status)

    if valid_flag && feasible_flag
        for j = 1:length(x._obbt_working_lower_index)
            if @inbounds x._obbt_working_lower_index[j]
                vi = @inbounds x._lower_variable_index[j]
                diff = MOI.get(x.relaxed_optimizer, MOI.VariablePrimal(), vi)
                diff -= @inbounds y.lower_variable_bounds[j]
                if abs(diff) <= x.obbt_tolerance
                    @inbounds x._obbt_working_lower_index[j] = false
                end
            end
        end
        for j = 1:length(x._obbt_working_upper_index)
            if @inbounds x._obbt_working_upper_index[j]
                vi = @inbounds x._lower_variable_index[j]
                diff = -MOI.get(x.relaxed_optimizer, MOI.VariablePrimal(), vi)
                diff += @inbounds y.upper_variable_bounds[j]
                if abs(diff) <= x.obbt_tolerance
                    @inbounds x._obbt_working_upper_index[j] = false
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
function aggressive_filtering!(x::Optimizer, y::NodeBB)

    # Initial filtering vector (negative one direction per remark in Gleixner2017)
    variable_number = x._variable_number
    v = -ones(variable_number)

    # Copy prior index set (ignores linear and binary terms)
    obbt_var_len = length(x.obbt_variable_values)
    copyto!(x._old_low_index, x._obbt_working_lower_index)
    copyto!(x._old_upp_index, x._obbt_working_upper_index)
    copyto!(x._new_low_index, x._obbt_working_lower_index)
    copyto!(x._new_upp_index, x._obbt_working_upper_index)

    # Exclude unbounded directions
    for i = 1:obbt_var_len
        if @inbounds x._new_low_index[i] && @inbounds y.lower_variable_bounds[i] == -Inf
            @inbounds x._new_low_index[i] = false
        end
    end
    for i = 1:obbt_var_len
        if @inbounds x._new_low_index[i] && @inbounds y.upper_variable_bounds[i] == Inf
            @inbounds x._new_low_index[i] = false
        end
    end

    # Begin the main algorithm
    for k = 1:x.obbt_aggressive_max_iteration

        # Set index differences and vector for filtering direction
        bool_indx_diff(x._lower_indx_diff, x._old_low_index, x._new_low_index)
        bool_indx_diff(x._upper_indx_diff, x._old_upp_index, x._new_upp_index)

        for i = 1:obbt_var_len
            if @inbounds x._lower_indx_diff[i] && @inbounds v[i] < 0.0
                @inbounds v[i] = 0.0
            end
        end
        for i = 1:obbt_var_len
            if @inbounds x._upper_indx_diff[i] && @inbounds v[i] > 0.0
                @inbounds v[i] = 0.0
            end
        end

        # Termination Condition
        ((~any(x._new_low_index) & ~any(x._new_upp_index)) || (iszero(v))) && break
        if k >= 2
            if (count(x._lower_indx_diff) + count(x._upper_indx_diff)) < x.obbt_aggressive_min_dimension
                break
            end
        end

        # Set objective in OBBT problem to filtering vector
        MOI.set(x.relaxed_optimizer, MOI.ObjectiveSense(), MOI.MAX_SENSE)
        saf = SAF(SAT.(v, x._lower_variable_index), 0.0)
        MOI.set(x.relaxed_optimizer, MOI.ObjectiveFunction{SAF}(), saf)

        # Optimizes the problem and if successful filter additional bounds
        MOI.optimize!(x.relaxed_optimizer)

        x._preprocess_termination_status = MOI.get(x.relaxed_optimizer, MOI.TerminationStatus())
        x._preprocess_result_status = MOI.get(x.relaxed_optimizer, MOI.PrimalStatus())
        valid_flag, feasible_flag = is_globally_optimal(x._preprocess_termination_status,
                                                        x._preprocess_result_status)

        if valid_flag
            if feasible_flag
                variable_primal = MOI.get(x.relaxed_optimizer, MOI.VariablePrimal(), x._lower_variable_index)
                copyto!(x._new_low_index, x._old_low_index)
                copyto!(x._new_upp_index, x._old_upp_index)
                for i = 1:obbt_var_len
                    if @inbounds x._old_low_index[i] && @inbounds variable_primal[i] == @inbounds y.lower_variable_bounds[i]
                        @inbounds x._new_low_index[i] = false
                    end
                end
                for i = 1:obbt_var_len
                    if @inbounds x._old_upp_index[i] && @inbounds variable_primal[i] == @inbounds y.upper_variable_bounds[i]
                        @inbounds x._new_upp_index[i] = false
                    end
                end
            end
        else
            return false
        end
    end
    copyto!(x._obbt_working_lower_index, x._new_low_index)
    copyto!(x._obbt_working_upper_index, x._new_upp_index)
    return true
end

"""
$(FUNCTIONNAME)

Performs OBBT with filtering and greedy ordering as detailed in:
Gleixner, A.M., Berthold, T., Müller, B. et al. J Glob Optim (2017) 67: 731.
https://doi.org/10.1007/s10898-016-0450-4
"""
function obbt(d::Optimizer)

    feasibility = true

    y = d._current_node

    println(" ")
    println("start obbt")
    println("y.lower_variable_bounds = $(y.lower_variable_bounds)")
    println("y.upper_variable_bounds = $(y.upper_variable_bounds)")

    @. d._current_xref = 0.5*(y.upper_variable_bounds + y.lower_variable_bounds)
    unsafe_check_fill!(isnan, d._current_xref, 0.0, length(d._current_xref))

    println(" ")
    println("d._current_xref = $(d._current_xref)")

    # solve initial problem to feasibility
    update_relaxed_problem_box!(d, y)
    relax_problem!(d, d._current_xref, 1)
    relax_objective!(d, d._current_xref)
    if d.objective_cut_on
        objective_cut_linear!(d, 1)
    end
    MOI.set(d.relaxed_optimizer, MOI.ObjectiveSense(), MOI.FEASIBILITY_SENSE)
    MOI.optimize!(d.relaxed_optimizer)

    # Sets indices to attempt OBBT on (full set...)
    obbt_variables = d.obbt_variable_values
    copyto!(d._obbt_working_lower_index, obbt_variables)
    copyto!(d._obbt_working_upper_index, obbt_variables)

    # Prefiltering steps && and sets initial LP values
    trivial_filtering!(d, y)
    if d.obbt_aggressive_on
        feasibility = aggressive_filtering!(d, y)
    end

    d._preprocess_termination_status = MOI.get(d.relaxed_optimizer, MOI.TerminationStatus())
    d._preprocess_result_status = MOI.get(d.relaxed_optimizer, MOI.PrimalStatus())
    valid_flag, feasible_flag = is_globally_optimal(d._preprocess_termination_status,
                                                    d._preprocess_result_status)

    if valid_flag && feasible_flag
        xLP = MOI.get(d.relaxed_optimizer, MOI.VariablePrimal(), d._lower_variable_index)
    else
        return false
    end

    while (any(d._obbt_working_lower_index) || any(d._obbt_working_upper_index)) & ~isempty(y)

        # Get lower value
        lower_indx = -1;     upper_indx = -1
        lower_value = Inf;   upper_value = Inf

        # min of xLP - yL on active
        if any(d._obbt_working_lower_index)
            for i = 1:length(d._obbt_working_lower_index)
                if @inbounds d._obbt_working_lower_index[i]
                    temp_value = @inbounds xLP[i] - y.lower_variable_bounds[i]
                    # Need less than or equal to handle unbounded cases
                    if temp_value <= lower_value
                        lower_value = temp_value
                        lower_indx = i
                    end
                end
            end
        end

        # min of yU - xLP on active
        if any(d._obbt_working_upper_index)
            for i = 1:length(d._obbt_working_upper_index)
                if @inbounds d._obbt_working_upper_index[i]
                    temp_value = @inbounds y.upper_variable_bounds[i] - xLP[i]
                    if temp_value <= upper_value
                        upper_value = temp_value
                        upper_indx = i
                    end
                end
            end
        end

        # default to upper bound if no lower bound is found, use maximum distance otherwise
        if lower_value <= upper_value && lower_indx > 0

            @inbounds d._obbt_working_lower_index[lower_indx] = false
            @inbounds var = d._lower_variable[lower_indx]
            MOI.set(d.relaxed_optimizer, MOI.ObjectiveSense(), MOI.MIN_SENSE)
            MOI.set(d.relaxed_optimizer, MOI.ObjectiveFunction{SV}(), var)
            MOI.optimize!(d.relaxed_optimizer)
            d._preprocess_termination_status = MOI.get(d.relaxed_optimizer, MOI.TerminationStatus())
            d._preprocess_result_status = MOI.get(d.relaxed_optimizer, MOI.PrimalStatus())
            valid_flag, feasible_flag = is_globally_optimal(d._preprocess_termination_status,
                                                            d._preprocess_result_status)

            if valid_flag
                if feasible_flag
                    xLP .= MOI.get(d.relaxed_optimizer, MOI.VariablePrimal(), d._lower_variable_index)
                    #(lower_indx === 1) && println("lower xLP[$lower_indx] = $(xLP[lower_indx])")
                    if is_integer_variable(d, lower_indx)
                        @inbounds y.lower_variable_bounds[lower_indx] = ceil(xLP[lower_indx])
                    else
                        @inbounds y.lower_variable_bounds[lower_indx] = xLP[lower_indx]
                    end
                    if isempty(y)
                        feasibility = false
                        break
                    end
                else
                    feasibility = false
                end
            else
                break
            end

        elseif upper_indx > 0

            d._obbt_working_upper_index[upper_indx] = false
            @inbounds var = d._lower_variable[upper_indx]
            MOI.set(d.relaxed_optimizer, MOI.ObjectiveSense(), MOI.MAX_SENSE)
            MOI.set(d.relaxed_optimizer, MOI.ObjectiveFunction{SV}(), var)
            MOI.optimize!(d.relaxed_optimizer)
            d._preprocess_termination_status = MOI.get(d.relaxed_optimizer, MOI.TerminationStatus())
            d._preprocess_result_status = MOI.get(d.relaxed_optimizer, MOI.PrimalStatus())
            valid_flag, feasible_flag = is_globally_optimal(d._preprocess_termination_status,
                                                            d._preprocess_result_status)

            if valid_flag
                if feasible_flag
                    xLP .= MOI.get(d.relaxed_optimizer, MOI.VariablePrimal(), d._lower_variable_index)
                    #(upper_indx === 1) && println("upper xLP[$upper_indx] = $(xLP[upper_indx])")
                    if is_integer_variable(d, upper_indx)
                        @inbounds y.upper_variable_bounds[upper_indx] = ceil(xLP[upper_indx])
                    else
                        @inbounds y.upper_variable_bounds[upper_indx] = min(xLP[upper_indx], y.upper_variable_bounds[upper_indx])
                    end
                    if isempty(y)
                        feasibility = false
                        break
                    end
                else
                    feasibility = false
                end
            else
                break
            end

        else
            break
        end
        trivial_filtering!(d, y)
    end
    println("end obbt")
    println("y.lower_variable_bounds = $(y.lower_variable_bounds)")
    println("y.upper_variable_bounds = $(y.upper_variable_bounds)")
    println(" ")

    return feasibility
end


"""

Performs feasibility-based bound tightening on a back-end constraint and returns `true` if it is feasible or
`false` if it is infeasible.
"""
function fbbt! end

function fbbt!(f::AffineFunctionIneq, n::NodeBB)

    # compute full sum
    lower_bounds = n.lower_variable_bounds
    upper_bounds = n.upper_variable_bounds
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
                (xh > xL) && (@inbounds lower_bounds[indx_k] = xh)
            elseif aik < 0.0
                (xh < xL) && return false
                (xh < xU) && (@inbounds upper_bounds[indx_k] = xh)
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

function fbbt!(f::AffineFunctionEq, n::NodeBB)
    # compute full sum
    lower_bounds = n.lower_variable_bounds
    upper_bounds = n.upper_variable_bounds
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
