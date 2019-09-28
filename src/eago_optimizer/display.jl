"""
    print_solution!

Prints solution information for the B&B problem. Displays first node found, solution value,
solution, and time spent solving subproblems.
"""
function print_solution!(x::Optimizer)
    if x.verbosity > 0
        println("First Solution Found at Node $(x._first_solution_node)")
        println("UBD = $(x._objective_value)")
        println("Solution is :")
        if (x._feasible_solution_found)
            for i=1:x._variable_number
                println("    X[$i] = $(x._continuous_solution[i])")
            end
        end
     end
     return
end

"""
    print_node!

Prints node information for the B&B problem. Node id, bound, and interval box.
"""
function print_node!(y::Optimizer)
    x = y._current_node
    println("Node ID: $(x.id), Lower Bound: $(x.lower_bound), Lower Variable Bounds:
             $(x.lower_variable_bounds), Upper Variable Bounds: $(x.upper_variable_bounds)")
    return
end

"""
    print_iteration!

Prints the iteration information based on verbosity. The header is displayed
every `header_interval`, the iteration info is displayed every `iteration_interval`.
"""
function print_iteration!(x::Optimizer)

    if (x.verbosity > 0)

        # prints header line every B.hdr_intv times
        if (mod(x._iteration_count, x.header_iterations) == 0 || x._iteration_count == 1)
            println("------------------------------------------------------------------------------------------------")
            println("|  Iteration #  |   Nodes  | Lower Bound | Upper Bound  |  Gap   |  Ratio  |  Time  | Time Left |")
            println("------------------------------------------------------------------------------------------------")
        end

        # prints iteration summary every B.itr_intv times
        if (mod(x._iteration_count, x.output_iterations) == 0)

            print_str = "| "

            max_len = 12
            temp_str = string(x._iteration_count)
            len_str = length(temp_str)
            print_str *= (" "^(max_len - len_str))*temp_str*"  | "

            max_len = 8
            temp_str = string(x._node_count)
            len_str = length(temp_str)
            print_str *= (" "^(max_len - len_str))*temp_str*" | "

            max_len = 11
            temp_str = string(round(x._global_lower_bound, sigdigits = 5))
            len_str = length(temp_str)
            print_str *= (" "^(max_len - len_str))*temp_str*" | "

            temp_str = string(round(x._global_upper_bound, sigdigits = 5))
            len_str = length(temp_str)
            print_str *= (" "^(max_len - len_str))*temp_str*"  |"

            max_len = 7
            temp_str = string(round(abs(x._global_upper_bound - x._global_lower_bound), sigdigits = 3))
            len_str = length(temp_str)
            print_str *= (" "^(max_len - len_str))*temp_str*" | "

            max_len = 7
            temp_str = string(round(relative_gap(x._global_lower_bound, x._global_upper_bound), sigdigits = 3))
            len_str = length(temp_str)
            print_str *= (" "^(max_len - len_str))*temp_str*" | "

            max_len = 7
            temp_str = string(round(x._run_time, sigdigits = 3))
            len_str = length(temp_str)
            print_str *= (" "^(max_len - len_str))*temp_str*" | "

            max_len = 7
            temp_str = string(round(x._time_left, sigdigits = 4))
            len_str = length(temp_str)
            print_str *= (" "^(max_len - len_str))*temp_str*"  |"

            println(print_str)
        end
    end

    return
end

"""
    print_results!

Prints the results of a single bounding problem.
"""
function print_results!(x::Optimizer, flag::Bool)
    if x.verbosity > 1
        if flag
            print("Lower Bound (First Iteration): $(x._lower_objective_value),")
            print(" Solution: $(x._lower_solution), Feasibility: $(x._lower_feasibility)\n")
        else
            print("Upper Bound: $(x._upper_objective_value), ")
            print(" Solution: $(x._upper_solution), Feasibility: $(x._upper_feasibility)\n")
        end
    end
    return
end

"""
    print_results_post_cut!

Prints the results after performing various cuts.
"""
function print_results_post_cut!(x::Optimizer)
    if x.verbosity > 1
        print("Lower Bound (Last Iteration): $(x._lower_objective_value)")
        print(", Solution: $(x._lower_solution), Feasibility: $(x._lower_feasibility)\n")
    end
    return
end
