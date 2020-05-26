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
# src/eago_optimizer/display.jl
# Functions used to print information about solution routine to console.
#############################################################################

"""
$(FUNCTIONNAME)

Prints solution information for the B&B problem. Displays first node found, solution value,
solution, and time spent solving subproblems.
"""
function print_solution!(m::Optimizer)
    if m._parameters.verbosity > 0
        println(" ")
        println("First Solution Found at Node $(m._first_solution_node)")
        println("UBD = $(MOI.get(m, MOI.ObjectiveValue()))")
        println("Solution is :")
        if m._feasible_solution_found
            for i = 1:m._input_problem._variable_number
                println("    X[$i] = $(m._continuous_solution[i])")
            end
        end
        println(" ")
     end
     return
end

"""
$(FUNCTIONNAME)

Prints node information for the B&B problem. Node id, bound, and interval box.
"""
function print_node!(m::Optimizer)
    x = m._current_node
    bound = (m._input_problem._optimization_sense === MOI.MIN_SENSE) ? x.lower_bound : -x.lower_bound
    println(" ")
    println("Node ID: $(x.id), Lower Bound: $(bound), Lower Variable Bounds:
             $(x.lower_variable_bounds), Upper Variable Bounds: $(x.upper_variable_bounds)")
    println(" ")
    return
end

const PRINTING_IOFORMAT = :SCI
const PRINTING_CHARSET = :ASCII

"""
$(FUNCTIONNAME)

Prints the iteration information based on verbosity. The header is displayed
every `header_interval`, the iteration info is displayed every `iteration_interval`.
"""
function print_iteration!(m::Optimizer)

    if m._parameters.verbosity > 0

        # prints header line every B.hdr_intv times
        if mod(m._iteration_count, m._parameters.header_iterations) === 0 || m._iteration_count === 1
            println("-----------------------------------------------------------------------------------------------------------------------------")
            println("|  Iteration #  |     Nodes    | Lower Bound  |  Upper Bound  |      Gap     |     Ratio    |     Time     |    Time Left   |")
            println("-----------------------------------------------------------------------------------------------------------------------------")
        end

        # prints iteration summary every B.itr_intv times
        if mod(m._iteration_count, m._parameters.output_iterations) === 0

            print_str = "| "

            max_len = 12
            temp_str = string(x._iteration_count)
            len_str = length(temp_str)
            print_str *= (" "^(max_len - len_str))*temp_str*"  | "

            max_len = 12
            temp_str = string(x._node_count)
            len_str = length(temp_str)
            print_str *= (" "^(max_len - len_str))*temp_str*" | "

            max_len = 12
            if x._input_problem._optimization_sense === MOI.MIN_SENSE
                lower = x._global_lower_bound
                upper = x._global_upper_bound
            else
                lower = -x._global_upper_bound
                upper = -x._global_lower_bound
            end
            #temp_str = string(round(lower, sigdigits = 5))
            #temp_str = string(lower, sigdigits = 3))
            temp_str = formatted(lower, PRINTING_IOFORMAT, ndigits=4, charset=PRINTING_CHARSET)
            len_str = length(temp_str)
            print_str *= (" "^(max_len - len_str))*temp_str*" | "

            #temp_str = formatted(upper, PRINTING_IOFORMAT, ndigits=4, charset=PRINTING_CHARSET)
            #temp_str = string(upper, sigdigits = 3))
            temp_str = formatted(upper, PRINTING_IOFORMAT, ndigits=4, charset=PRINTING_CHARSET)
            len_str = length(temp_str)
            print_str *= (" "^(max_len - len_str))*temp_str*"  |"

            max_len = 12
            #temp_str = string(round(abs(x._global_upper_bound - x._global_lower_bound), sigdigits = 3))
            temp_str = formatted(abs(x._global_upper_bound - x._global_lower_bound), PRINTING_IOFORMAT, ndigits=4, charset=PRINTING_CHARSET)
            len_str = length(temp_str)
            print_str *= (" "^(max_len - len_str))*temp_str*"  | "

            max_len = 12
            #temp_str = string(round(relative_gap(x._global_lower_bound, x._global_upper_bound), sigdigits = 3))
            temp_str = formatted(relative_gap(x._global_lower_bound, x._global_upper_bound), PRINTING_IOFORMAT, ndigits=4, charset=PRINTING_CHARSET)
            len_str = length(temp_str)
            print_str *= (" "^(max_len - len_str))*temp_str*" | "

            max_len = 12
            #temp_str = string(round(x._run_time, sigdigits = 3))
            temp_str = formatted(x._run_time, PRINTING_IOFORMAT, ndigits=4, charset=PRINTING_CHARSET)
            len_str = length(temp_str)
            print_str *= (" "^(max_len - len_str))*temp_str*" |  "

            max_len = 12
            #temp_str = string(round(x._time_left, sigdigits = 4))
            temp_str = formatted(x._time_left, PRINTING_IOFORMAT, ndigits=4, charset=PRINTING_CHARSET)
            len_str = length(temp_str)
            print_str *= (" "^(max_len - len_str))*temp_str*"  |"

            println(print_str)
        end
    end

    return
end

"""
$(FUNCTIONNAME)

Prints the results of a single bounding problem.
"""
function print_results!(m::Optimizer, flag::Bool)
    if m._parameters.verbosity > 1
        println(" ")
        if flag
            obj_val = m._lower_objective_value
            if m._input_problem._optimization_sense === MOI.MIN_SENSE
                print("Lower Bound (First Iteration): $(obj_val),")
            else
                print("Upper Bound (First Iteration): $(-obj_val),")
            end
            print(" Solution: $(m._lower_solution), Feasibility: $(m._lower_feasibility)\n")
            println("Termination Status Code: $(m._lower_termination_status)")
            println("Result Code: $(m._lower_result_status)")
        else
            obj_val = m._upper_objective_value
            if m._input_problem._optimization_sense === MOI.MIN_SENSE
                print("Upper Bound: $(obj_val), ")
            else
                print("Lower Bound: $(-obj_val), ")
            end
            print(" Solution: $(m._upper_solution), Feasibility: $(m._upper_feasibility)\n")
            println("Termination Status Code: $(m._upper_termination_status)")
            println("Result Code: $(m._upper_result_status)")
        end
        println(" ")
    end
    return
end

"""
$(FUNCTIONNAME)

Prints the results after performing various cuts.
"""
function print_results_post_cut!(m::Optimizer)
    if x._parameters.verbosity > 1
        println(" ")
        if m._input_problem._optimization_sense === MOI.MIN_SENSE
            print("Lower Bound (Last Iteration): $(m._lower_objective_value)")
        else
            print("Upper Bound (Last Iteration): $(-m._lower_objective_value)")
        end
        print(", Solution: $(m._lower_solution), Feasibility: $(m._lower_feasibility)\n")
        println(" ")
    end
    return
end
