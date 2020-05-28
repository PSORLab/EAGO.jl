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
# Printing is done with reference to the input problem is there is any
# ambiguity.
#############################################################################

const PRINTING_IOFORMAT = :SCI
const PRINTING_CHARSET = :ASCII

"""
$(FUNCTIONNAME)
"""
function display_relaxed_optimizer!(m::Optimizer, optimizer::T, note::String) where T

    println("------------------------------------------------------------------------------------------")
    println("                                           ")
    println("Displaying optimizer"*note)
    println("                                           ")

    objective_function_type = MOI.get(optimizer, MOI.ObjectiveFunctionType())
    if objective_function_type == SV
        objective_function = MOI.get(optimizer, MOI.ObjectiveFunction{SV}())
        variable_index_value = objective_function.variable.value
        println("Objective function (Single Variable) = x[$(variable_index_value)]")
    elseif objective_function_type == SAF
        objective_function = MOI.get(optimizer, MOI.ObjectiveFunction{SAF}())
        println("Objective function (Scalar Affine) = x[$(objective_function)]")
    end

    objective_sense = MOI.get(optimizer, MOI.ObjectiveSense())
    println("Objective sense isa $(objective_sense) \n")

    variable_number = MOI.get(optimizer, MOI.ListOfVariableIndices())
    variable_info = Vector{Any}[Any[nothing for j = 1:5] for i in variable_number]
    sv_lt = MOI.get(optimizer, MOI.ListOfConstraintIndices{SV, LT}())
    for temp1 in sv_lt
        sv_lt_func = MOI.get(optimizer, MOI.ConstraintFunction(), temp1)
        sv_lt_set = MOI.get(optimizer, MOI.ConstraintSet(), temp1)
        variable_info[sv_lt_func.variable.value][1] = sv_lt_func.variable.value
        variable_info[sv_lt_func.variable.value][2] = "LessThan    "
        variable_info[sv_lt_func.variable.value][5] = sv_lt_set.upper
    end

    sv_gt = MOI.get(optimizer, MOI.ListOfConstraintIndices{SV, GT}())
    for temp2 in sv_gt
        sv_gt_func = MOI.get(optimizer, MOI.ConstraintFunction(), temp2)
        sv_gt_set = MOI.get(optimizer, MOI.ConstraintSet(), temp2)
        variable_info[sv_gt_func.variable.value][1] = sv_gt_func.variable.value
        if variable_info[sv_gt_func.variable.value][2] == "LessThan    "
            variable_info[sv_gt_func.variable.value][2] = "Interval    "
        else
            variable_info[sv_gt_func.variable.value][2] = "GreaterThan "
        end
        variable_info[sv_gt_func.variable.value][4] = sv_gt_set.lower
    end

    sv_et = MOI.get(optimizer, MOI.ListOfConstraintIndices{SV, ET}())
    for temp3 in sv_et
        sv_et_func = MOI.get(optimizer, MOI.ConstraintFunction(), temp3)
        sv_et_set = MOI.get(optimizer, MOI.ConstraintSet(), temp3)
        variable_info[sv_et_func.variable.value][1] = sv_et_func.variable.value
        variable_info[sv_et_func.variable.value][2] = "EqualTo     "
        variable_info[sv_et_func.variable.value][4] = sv_et_set.value
        variable_info[sv_et_func.variable.value][5] = sv_et_set.value
    end

    sv_zo = MOI.get(optimizer, MOI.ListOfConstraintIndices{SV, ZO}())
    for temp4 in sv_zo
        sv_zo_func = MOI.get(optimizer, MOI.ConstraintFunction(), temp4)
        sv_zo_set = MOI.get(optimizer, MOI.ConstraintSet(), tem4)
        variable_info[sv_et_func.variable.value][1] = sv_zo_func.variable.value
        variable_info[sv_et_func.variable.value][2] = "ZeroOne     "
    end

    count = 1
    for i = 1:length(variable_info)
        if variable_info[i][1] == nothing
            variable_info[i][1] = i
            variable_info[i][2] = "            "
            variable_info[i][4] = "            "
            variable_info[i][5] = "            "
        end
        if m._branch_variables[i]
            nindx = m._sol_to_branch_map[i]
            variable_info[i][3] = m._current_xref[nindx]
        else
            variable_info[i][3] = "             "
        end
    end

    println("------------------------------------------------------------------------------------------")
    println("| Table of Variables                                                                     |")
    println("------------------------------------------------------------------------------------------")
    println("|    Variable   |     Bounds Type    |    Last Cut    |   Lower Bound   |   Upper Bound   |")
    println("------------------------------------------------------------------------------------------")

    for (i, var) in enumerate(variable_number)

        print_str = "|     "
        (i < 10) && (print_str = print_str*" ")
        (i < 100) && (print_str = print_str*" ")
        (i < 1000) && (print_str = print_str*" ")

        print_str = print_str*"x[$(variable_info[i][1])]   |"   # Index number
        print_str = print_str*"       $(variable_info[i][2]) |"  # Bounds type

        if variable_info[i][3] isa Number
            print_str = print_str*"       "*formatted(variable_info[i][3], PRINTING_IOFORMAT, ndigits=4, charset=PRINTING_CHARSET)*" |"    # Lower bound if any
        else
            print_str = print_str*"                |"
        end

        if variable_info[i][4] isa Number
            print_str = print_str*"       "*formatted(variable_info[i][4], PRINTING_IOFORMAT, ndigits=4, charset=PRINTING_CHARSET)*" |"    # Lower bound if any
        else
            print_str = print_str*"               |"
        end

        if variable_info[i][5] isa Number
            print_str = print_str*"       "*formatted(variable_info[i][5], PRINTING_IOFORMAT, ndigits=4, charset=PRINTING_CHARSET)*"    |"    # Upper bound if any
        else
            print_str = print_str*"                  | "
        end

        println(print_str)

    end
    println("------------------------------------------------------------------------------------------\n")

    println("------------------------------------------------------------------------------------------")
    println("| Table of Affine Relaxations                                                            |")
    println("------------------------------------------------------------------------------------------")
    header_str = "| source type   | source loc   | affine relaxation     |"
    println(header_str)

    saf_lt = MOI.get(optimizer, MOI.ListOfConstraintIndices{SAF, LT}())
    for temp5 in saf_lt
        saf = " "
        saf_lt_func = MOI.get(optimizer, MOI.ConstraintFunction(), temp5)
        for temp6 in saf_lt_func.terms
            saf = saf*"$(temp6.coefficient)*x[$(temp6.variable_index.value)] +"
        end
        saf_lt_set = MOI.get(optimizer, MOI.ConstraintSet(), temp5)
        saf = saf*formatted(saf_lt_set.upper, PRINTING_IOFORMAT, ndigits=4, charset=PRINTING_CHARSET)
        saf = saf*" <= 0.0"
        println(saf)
    end


    println("------------------------------------------------------------------------------------------")

    #println("constant = $(saf.constant)")
    #str = "["
    #for term in saf.terms
    #    coeff = term.coefficient
    #    index = term.variable_index
    #    str *= " ($coeff, $index) "
    #end
    #str *= " ]"
    #println(str)
    println("------------------------------------------------------------------------------------------")
    nothing
end

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
            for i = 1:m._input_problem._variable_count
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
    n = m._current_node
    bound = (m._input_problem._optimization_sense === MOI.MIN_SENSE) ? n.lower_bound : -n.lower_bound
    println(" ")
    println("Node ID: $(n.id), Lower Bound: $(bound), Lower Variable Bounds:
             $(n.lower_variable_bounds), Upper Variable Bounds: $(n.upper_variable_bounds)")
    println(" ")
    return
end

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
            temp_str = string(m._iteration_count)
            len_str = length(temp_str)
            print_str *= (" "^(max_len - len_str))*temp_str*"  | "

            max_len = 12
            temp_str = string(m._node_count)
            len_str = length(temp_str)
            print_str *= (" "^(max_len - len_str))*temp_str*" | "

            max_len = 12
            if m._input_problem._optimization_sense === MOI.MIN_SENSE
                lower = m._global_lower_bound
                upper = m._global_upper_bound
            else
                lower = -m._global_upper_bound
                upper = -m._global_lower_bound
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
            temp_str = formatted(abs(m._global_upper_bound - m._global_lower_bound), PRINTING_IOFORMAT, ndigits=4, charset=PRINTING_CHARSET)
            len_str = length(temp_str)
            print_str *= (" "^(max_len - len_str))*temp_str*"  | "

            max_len = 12
            #temp_str = string(round(relative_gap(x._global_lower_bound, x._global_upper_bound), sigdigits = 3))
            temp_str = formatted(relative_gap(m._global_lower_bound, m._global_upper_bound), PRINTING_IOFORMAT, ndigits=4, charset=PRINTING_CHARSET)
            len_str = length(temp_str)
            print_str *= (" "^(max_len - len_str))*temp_str*" | "

            max_len = 12
            #temp_str = string(round(x._run_time, sigdigits = 3))
            temp_str = formatted(m._run_time, PRINTING_IOFORMAT, ndigits=4, charset=PRINTING_CHARSET)
            len_str = length(temp_str)
            print_str *= (" "^(max_len - len_str))*temp_str*" |  "

            max_len = 12
            #temp_str = string(round(x._time_left, sigdigits = 4))
            temp_str = formatted(m._time_left, PRINTING_IOFORMAT, ndigits=4, charset=PRINTING_CHARSET)
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
    if m._parameters.verbosity > 1
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
