"""
    solve_nlp!

Solves the branch and bound problem with the input EAGO optimizer object.
"""
function global_solve!(x::Optimizer)

  # initializes Flags
  if (~x.warm_start)
    x.current_iteration_count = 1
    x.current_node_count = 0
  end

  PreprocessTime = 0.0
  LowerProblemTime = 0.0
  UpperProblemTime = 0.0
  PostprocessTime = 0.0

  ran_preprocess = false
  ran_lower_problem = false
  ran_upper_problem = false
  ran_postprocess = false
  conv_check = false

  # terminates when max nodes or iteration is reach, or when node stack is empty
  iterationcountinternal = 0
  while termination_check(x)

    iterationcountinternal += 1

    LowerProblemTime = 0.0
    UpperProblemTime = 0.0
    PostprocessTime = 0.0

    ran_lower_problem = false
    ran_upper_problem = false
    conv_check = false

    x._cut_iterations = 0

    # Fathom nodes with lower bound greater than global upper bound
    x.global_lower_bound = find_lower_bound(x)
    x.history.lower_bound[x.current_iteration_count] = x.global_lower_bound

     # Selects node, deletes it from stack, prints based on verbosity
    CurrentKey,CurrentNode = node_selection(x);  x.current_node_count -= 1
    (x.verbosity >= 3) && print_node!(CurrentKey,CurrentNode) # Prints node in full verbosity mode

    # Solves preprocessing/LBD/UBD/postprocessing once to get timing right
    x.current_preprocess_info.feasibility = true
    x.current_postprocess_info.feasibility = true

    #=
    if (x.current_iteration_count == 1)
      OldLowerInfo = LowerInfo(x.current_lower_info)
      OldUpperInfo = UpperInfo(x.current_upper_info)
      OldPreprocessInfo = PreprocessInfo(x.current_preprocess_info)
      OldPostprocessInfo = PostprocessInfo(x.current_postprocess_info)
      tempNode = NodeBB(CurrentNode)
      (x.verbosity >= 4) && println("started initial preprocessing")
      preprocess!(x,tempNode)
      (x.verbosity >= 4) && println("finished initial preprocessing")
      lower_problem!(x,tempNode)
      (x.verbosity >= 4) && println("finished initial lower problem")
      upper_problem!(x,tempNode)
      (x.verbosity >= 4) && println("finished initial upper problem")
      postprocess!(x,tempNode)
      (x.verbosity >= 4) && println("finished initial postprocessing")
      x.current_lower_info = OldLowerInfo
      x.current_upper_info = OldUpperInfo
      x.current_preprocess_info = OldPreprocessInfo
      x.current_postprocess_info = OldPostprocessInfo
    end
    x.current_preprocess_info.feasibility = true
    x.current_postprocess_info.feasibility = true
    =#

    # Performs prepocessing and times
    PreprocessTime = @elapsed preprocess!(x,CurrentNode)

    x.current_upper_info.feasibility = true

    ran_lower_problem = x.current_preprocess_info.feasibility
    if ran_lower_problem
      # solves & times lower bounding problem
      LowerProblemTime = @elapsed lower_problem!(x,CurrentNode)
      x.history.lower_bound[x.current_iteration_count] = x.global_lower_bound
      print_results!(x,true)

      while cut_condition(x)
        LowerProblemTime += @elapsed add_cut!(x, CurrentNode)
        x._cut_iterations += 1
        x.history.lower_bound[x.current_iteration_count] = x.global_lower_bound
      end
      print_results_post_cut!(x)
      x.history.cut_count[x.current_iteration_count] = x._cut_iterations

      # checks for infeasibility stores solution
      ran_upper_problem = x.current_lower_info.feasibility
      if ran_upper_problem
        conv_check = ~convergence_check(x)
        if conv_check
          UpperProblemTime = @elapsed upper_problem!(x,CurrentNode)
          print_results!(x,false)

          # Stores information if better feasible upper bound is formed
          ran_postprocess = x.current_upper_info.feasibility
          if ran_postprocess
            if (x.current_upper_info.value < x.global_upper_bound)
              x.feasible_solution_found = true
              x.first_solution_node = x.maximum_node_id
              x.solution_value = x.current_upper_info.value
              x.continuous_solution[:] = x.current_upper_info.solution
              x.history.upper_bound[x.current_iteration_count] = x.solution_value
              if (x.optimization_sense == MOI.FEASIBILITY_SENSE) && (~x.feasible_local_continue || x.local_solve_only)
                break
              end
            else
              x.history.upper_bound[x.current_iteration_count] = x.history.upper_bound[x.current_iteration_count-1]
            end
          else
            x.history.upper_bound[x.current_iteration_count] = x.history.upper_bound[x.current_iteration_count-1]
          end

          # Performs and times post processing
          PostprocessTime = @elapsed postprocess!(x,CurrentNode)

          # Checks to see if the node
          if (x.current_postprocess_info.feasibility)
            if repeat_check(x, CurrentNode)
              single_storage!(x, CurrentNode)
              x.node_repetitions += 1
            else
              Y1,Y2 = bisection_function(x,CurrentNode)
              node_storage!(x,Y1,Y2)
              x.node_repetitions = 1
            end
          end
        end
      end
      fathom!(x)
    else
      x.current_lower_info.value = -Inf
      x.current_lower_info.feasibility = false
      x.current_upper_info.feasibility = false
    end

    ~ran_lower_problem && (x.history.lower_bound[x.current_iteration_count] = x.history.lower_bound[x.current_iteration_count-1])
    (~ran_upper_problem) && (x.history.upper_bound[x.current_iteration_count] = x.history.upper_bound[x.current_iteration_count-1])
    (ran_upper_problem && ~conv_check) && (x.history.upper_bound[x.current_iteration_count] = x.history.upper_bound[x.current_iteration_count-1])

    x.history.preprocess_time[x.current_iteration_count] = x.history.preprocess_time[x.current_iteration_count-1] + PreprocessTime
    x.history.lower_time[x.current_iteration_count] = x.history.lower_time[x.current_iteration_count-1] + LowerProblemTime
    x.history.upper_time[x.current_iteration_count] = x.history.upper_time[x.current_iteration_count-1] + UpperProblemTime
    x.history.postprocess_time[x.current_iteration_count] = x.history.postprocess_time[x.current_iteration_count-1] + PostprocessTime

    x.history.count[x.current_iteration_count] = x.current_node_count
    x.history.lower_count += ran_lower_problem
    x.history.upper_count += ran_upper_problem

    print_iteration!(x)
    x.current_iteration_count += 1
  end

  x.objective_value = x.global_upper_bound
  if (x.optimization_sense == MOI.MAX_SENSE)
        x.objective_value *= -1.0
  end
  print_solution!(x)                              # Prints the solution
end
