"""
    print_solution!

Prints solution information for the B&B problem. Displays first node found, solution value,
solution, and time spent solving subproblems.
"""
function print_solution!(x::Optimizer)
  temp::Float64 = 0.0
  if (x.verbosity > 0)
    println("First Solution Found at Node $(x.first_solution_node)")
    println("UBD = $(x.objective_value)")
    println("Solution is :")
    if (x.feasible_solution_found)
      xlen = length(x.continuous_solution)
      xlen += (x.reform_epigraph_flag) ? -1 : 0
      for i=1:xlen
        temp = x.continuous_solution[i]
        println("    X[$i] = $temp")
      end
    end
    #=
    println("Total LBD problems solved = $(x.History.LowerCount) in $(x.History.LowerTime[x.CurrentIterationCount]) seconds.")
    println("Total UBD problems solved = $(x.History.UpperCount) in $(x.History.UpperTime[x.CurrentIterationCount]) seconds.")
    println("Total time spent preprocessing =  $(x.History.PreprocessTime[x.CurrentIterationCount]) seconds.")
    println("Total time spent postprocessing = $(x.History.PostprocessTime[x.CurrentIterationCount]) seconds.")
    =#
  end
end

"""
    print_node!

Prints node information for the B&B problem. Node id, bound, and interval box.
"""
function print_node!(id::Int,x::NodeBB)
    println("Node ID: $(id), Lower Bound: $(x.lower_bound), Lower Variable Bounds:
             $(x.lower_variable_bounds), Upper Variable Bounds: $(x.upper_variable_bounds)")
end


"""
    print_iteration!

Prints the iteration information based on verbosity. The header is displayed
every `header_interval`, the iteration info is displayed every `iteration_interval`.
"""
function print_iteration!(x::Optimizer)
  if (x.verbosity > 0)
    # prints header line every B.hdr_intv times
    if (mod(x.current_iteration_count,x.header_iterations) == 0 || x.current_iteration_count == 1)
      println("Iteration   NodeID    Current_LBD     Global_LBD     Global_UBD      NodesLeft     Absolute_Gap    Absolute_Ratio     LBD_Feas     UBD_Feas")
    end
    # prints iteration summary every B.itr_intv times
    sbool1 = x.current_lower_info.feasibility ? "true" : "false"
    sbool2 = x.current_upper_info.feasibility ? "true" : "false"
    nid = x.global_upper_bound
    lbdp = x.current_lower_info.feasibility ? x.current_lower_info.value : Inf
    if ((mod(x.current_iteration_count,x.output_iterations) == 0))
      ptr_arr1 = join([@sprintf("%6u",x) for x in Int[x.current_iteration_count x.current_node_count]], ",   ")
      ptr_arr2 = join([@sprintf("%3.7f",x) for x in Float64[lbdp x.global_lower_bound x.global_upper_bound]], ",     ")
      ptr_arr3 = join([@sprintf("%6u",x) for x in Int[x.current_node_count]], ",")
      ptr_arr4 = join([@sprintf("%3.7f",x) for x in Float64[abs(x.global_upper_bound-x.global_lower_bound),
                            abs(x.global_upper_bound-x.global_lower_bound)/(min(abs(x.global_lower_bound),abs(x.global_upper_bound)))]], ",       ")
      ptr_arr5 = join([@sprintf("%s",x) for x in String[sbool1 sbool2]], ",       ")
#      ptr_arr1 = join([@sprintf("%6u",x) for x in Int[k_int nid]], ",   ")
#      ptr_arr2 = join([@sprintf("%3.7f",x) for x in Float64[lbdp lbd ubd]], ",     ")
#      ptr_arr3 = join([@sprintf("%6u",x) for x in Int[k_nod]], ",")
#      ptr_arr4 = join([@sprintf("%3.7f",x) for x in Float64[abs(ubd-lbd) abs(ubd-lbd)/(min(abs(lbd),abs(ubd)))]], ",       ")
#      ptr_arr5 = join([@sprintf("%s",x) for x in Bool[sbool1 sbool2]], ",       ")
      println(string(ptr_arr1,",      ",ptr_arr2,",      ",ptr_arr3),",        ", ptr_arr4,",        ", ptr_arr5)
    end
  end
end

"""
    print_results!

Prints the results of a single bounding problem.
"""
function print_results!(B::Optimizer,lbd_bool::Bool)
  if (B.verbosity > 1)
    if (lbd_bool)
      println("Lower Bound (First Iteration): $(B.current_lower_info.value), Solution: $(B.current_lower_info.solution), Feasibility: $(B.current_lower_info.feasibility)")
    else
      println("Upper Bound: $(B.current_upper_info.value), Solution: $(B.current_upper_info.solution), Feasibility: $(B.current_upper_info.feasibility)")
    end
  end
end

function print_results_post_cut!(B::Optimizer)
  if (B.verbosity > 1)
    println("Lower Bound (Last Iteration): $(B.current_lower_info.value), Solution: $(B.current_lower_info.solution), Feasibility: $(B.current_lower_info.feasibility)")
  end
end
