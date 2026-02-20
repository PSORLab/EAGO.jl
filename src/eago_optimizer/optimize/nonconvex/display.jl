# Copyright (c) 2018 Matthew Wilhelm, Robert Gottlieb, Dimitri Alston,
# Alireza Miraliakbar, Matthew Stuber, and the University of Connecticut (UConn)
# This code is licensed under the MIT license (see LICENSE.md for full details).
################################################################################
# EAGO
# A development environment for robust and global optimization
# https://github.com/PSORLab/EAGO.jl
################################################################################
# src/eago_optimizer/optimize/nonconvex/display.jl
# Functions used to print information about solution routine to console.
# Printing is done with reference to the input problem if there is any
# ambiguity.
################################################################################

"""
$(FUNCTIONNAME)

Print solution information for the B&B problem. Display node with the best solution,
solution value, solution, and time spent solving subproblems. This print occurs
following termination of the B&B algorithm.
"""
function print_solution!(m::GlobalOptimizer)
    if _verbosity(m) > 0
        println(" ")
        if m._end_state == GS_OPTIMAL
            println("Empty Stack: Exhaustive Search Finished")
        elseif m._end_state == GS_INFEASIBLE
            println("Empty Stack: Infeasible")
        elseif m._end_state == GS_NODE_LIMIT
            println("Node Limit Exceeded")
        elseif m._end_state == GS_ITERATION_LIMIT
            println("Iteration Limit Exceeded")
        elseif m._end_state == GS_RELATIVE_TOL
            println("Relative Tolerance Achieved")
        elseif m._end_state == GS_ABSOLUTE_TOL
            println("Absolute Tolerance Achieved")
        elseif m._end_state == GS_TIME_LIMIT
            println("Time Limit Exceeded")
        end
        if m._end_state == GS_OPTIMAL || m._end_state == GS_RELATIVE_TOL || m._end_state == GS_ABSOLUTE_TOL
            println("Optimal Solution Found at Node $(m._solution_node)")
            println("Lower Bound: $(MOI.get(m, MOI.ObjectiveBound()))")
            println("Upper Bound: $(MOI.get(m, MOI.ObjectiveValue()))")
        elseif m._end_state == GS_INFEASIBLE
            println("No Solution Found")
        else
            println("Best Solution Found at Node $(m._solution_node)")
            println("Lower Bound: $(MOI.get(m, MOI.ObjectiveBound()))")
            println("Upper Bound: $(MOI.get(m, MOI.ObjectiveValue()))")
        end
        if m._feasible_solution_found
            println("Solution:")
            variable_names = String[]
            for i = 1:m._input_problem._variable_count
                push!(variable_names, string(m._input_problem._variable_names[MOI.VariableIndex(i)]))
            end
            maxlen = maximum(length.(variable_names))
            addlen = maxlen .- length.(variable_names)
            print_list = " ".^addlen.*variable_names
            for i = 1:m._input_problem._variable_count
                println("   $(print_list[i]) = $(m._continuous_solution[i])")
            end
        end
        println(" ")
     end
     return
end

"""
$(FUNCTIONNAME)

Print information about the current node. Includes node ID, lower bound,
upper bound, and interval box.
"""
function print_node!(m::GlobalOptimizer)
    if _verbosity(m) >= 3
        n = m._current_node
        lower_bound = _is_input_min(m) ? n.lower_bound : -n.lower_bound
        upper_bound = _is_input_min(m) ? n.upper_bound : -n.upper_bound
        k = length(n) - (_obj_var_slack_added(m) ? 1 : 0)
        println(" ")
        println("Node ID: $(n.id)")
        println("Lower Bound: $(lower_bound)")
        println("Upper Bound: $(upper_bound)")
        println("Lower Variable Bounds: $(n.lower_variable_bounds[1:k])")
        println("Upper Variable Bounds: $(n.upper_variable_bounds[1:k])")
        println(" ")
    end
    return
end

"""
$(FUNCTIONNAME)

Print status information based on iteration count. The header print frequency is
based on the `header_iterations` setting, and the data print frequency is based on
the `output_iterations` setting.
"""
function print_iteration!(m::GlobalOptimizer, end_flag::Bool)

    if _verbosity(m) > 0

        # Print header line every `header_iterations` times and print iteration summary every `output_iterations` times
        if m._last_printed_iteration != m._iteration_count && (mod(m._iteration_count, m._parameters.output_iterations) === 0 || end_flag)
            if m._iteration_count == m._parameters.output_iterations || mod(m._iteration_count, m._parameters.header_iterations) < m._parameters.output_iterations
                println("-----------------------------------------------------------------------------------------------------------------")
                println("| Iteration # |    Nodes    | Lower Bound | Upper Bound |     Gap     |    Ratio    |    Timer    |  Time Left  |")
                println("-----------------------------------------------------------------------------------------------------------------")
            end
            # Print start
            print_str = "| "

            # Print iteration number
            max_len = 11
            temp_str = string(m._iteration_count)
            len_str = length(temp_str)
            print_str *= (" "^(max_len - len_str))*temp_str*" | "

            # Print node count
            max_len = 11
            temp_str = string(m._node_count)
            len_str = length(temp_str)
            print_str *= (" "^(max_len - len_str))*temp_str*" | "

            # Print lower bound
            max_len = 11
            temp_str = @sprintf "%.3E" m._global_lower_bound
            len_str = length(temp_str)
            print_str *= (" "^(max_len - len_str))*temp_str*" | "

            # Print upper bound
            max_len = 11
            temp_str = @sprintf "%.3E" m._global_upper_bound
            len_str = length(temp_str)
            print_str *= (" "^(max_len - len_str))*temp_str*" | "

            # Print absolute gap between lower and upper bound
            max_len = 11
            temp_str = @sprintf "%.3E" abs(m._global_upper_bound - m._global_lower_bound)
            len_str = length(temp_str)
            print_str *= (" "^(max_len - len_str))*temp_str*" | "

            # Print relative gap between lower and upper bound
            max_len = 11
            temp_str = @sprintf "%.3E" relative_gap(m._global_lower_bound, m._global_upper_bound)
            len_str = length(temp_str)
            print_str *= (" "^(max_len - len_str))*temp_str*" | "

            # Print run time
            max_len = 11
            temp_str = @sprintf "%.2F" m._run_time
            len_str = length(temp_str)
            print_str *= (" "^(max_len - len_str))*temp_str*" | "

            # Print time remaining
            max_len = 11
            temp_str = @sprintf "%.2F" m._time_left
            len_str = length(temp_str)
            print_str *= (" "^(max_len - len_str))*temp_str*" |"

            println(print_str)

            # Update printed iteration
            m._last_printed_iteration = m._iteration_count
        end
        if end_flag
            println("-----------------------------------------------------------------------------------------------------------------")
        end
    end

    return
end

"""
$(TYPEDSIGNATURES)

Print the results of a single (lower or upper) bounding problem. `lower_flag=true`
prints information for the lower problem, `lower_flag=false` prints information for
the upper problem.
"""
function print_results!(m::GlobalOptimizer, lower_flag::Bool)
    if _verbosity(m) >= 2
        k = length(m._lower_solution) - (_obj_var_slack_added(m) ? 1 : 0)
        println(" ")
        if lower_flag
            if _is_input_min(m)
                println("Lower Bound (First Iteration): $(m._lower_objective_value)")
            else
                println("Upper Bound (First Iteration): $(m._lower_objective_value)")
            end
            println("Solution: $(m._lower_solution[1:k])")
            println("Feasibility: $(m._lower_feasibility)")
            println("Termination Status Code: $(m._lower_termination_status)")
            println("Result Code: $(m._lower_primal_status)")
        else
            if _is_input_min(m)
                println("Upper Bound: $(m._upper_objective_value)")
            else
                println("Lower Bound: $(m._upper_objective_value)")
            end
            println("Solution: $(m._upper_solution[1:k])")
            println("Feasibility: $(m._upper_feasibility)")
            println("Termination Status Code: $(m._upper_termination_status)")
            println("Result Code: $(m._upper_result_status)")
        end
        println(" ")
    end
    return
end

"""
$(FUNCTIONNAME)

Print noteworthy information prior to running branch-and-bound. Currently prints
a note about flipping `max(f)` to `-min(-f)` internally, if the input is a 
maximization problem and `verbosity >= 3`.
"""
function print_preamble!(m::GlobalOptimizer)
    if _verbosity(m) >= 3
        if !_is_input_min(m) && iszero(m._iteration_count)
            @info("""
            For maximization problems, a transformation from `max(f)` to `-min(-f)` is applied.
            Objective values for each subproblem are the negative value of the objective
            in the original problem and are reconciled after branch-and-bound terminates.""")
        end
    end
    return
end
