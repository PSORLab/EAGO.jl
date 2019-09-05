"""
    default_cut_condition

Branch-and-cut feature currently under development. Currently, returns false.
"""
cut_condition(x::Optimizer) = x._cut_add_flag && (x._cut_iterations < x.cut_max_iterations)

#=
function check_cut_tolerance(x::Optimizer, solution::Vector{Float64})
end
=#

"""
    default_add_cut!

Branch-and-Cut under development.
"""
function add_cut!(x::Optimizer, y::NodeBB)

    xpnt = x.current_lower_info.solution[1:end]
    update_lower_variable_bounds1!(x, y, x.working_relaxed_optimizer)
    x.relax_function!(x, x.working_relaxed_optimizer, y, x.relaxation, xpnt, load = true)
    x.relax_function!(x, x.working_relaxed_optimizer, y, x.relaxation, xpnt, load = false)

    # Optimizes the object
    MOI.optimize!(x.working_relaxed_optimizer)

    # Process output info and save to CurrentUpperInfo object
    termination_status = MOI.get(x.working_relaxed_optimizer, MOI.TerminationStatus())
    result_status_code = MOI.get(x.working_relaxed_optimizer, MOI.PrimalStatus())
    valid_flag, feasible_flag = is_globally_optimal(termination_status, result_status_code)
    last_obj = x.current_lower_info.value
    #println("cut solution: $last_obj")
    if valid_flag
        if feasible_flag
            x.current_lower_info.feasibility = true
            objval = MOI.get(x.working_relaxed_optimizer, MOI.ObjectiveValue())
            x.current_lower_info.value = objval
            #interval_lower_bound!(x, y, false)
            vprimal_solution = MOI.get(x.working_relaxed_optimizer, MOI.VariablePrimal(), x.lower_variables)
            solutions_distinct = ~isapprox(last_obj, objval, atol = x.absolute_tolerance/2.0)
            x.cut_add_flag = solutions_distinct
            x.current_lower_info.solution[1:end] = vprimal_solution
            set_dual!(x)
        else
            x.cut_add_flag = false
            x.current_lower_info.feasibility = false
            x.current_lower_info.value = -Inf
        end
    else
        x.cut_add_flag = false
    end

end
