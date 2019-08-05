function local_nlp_solve!(x::Optimizer)
  current_iteration_count = 1
  CurrentKey, CurrentNode = x.node_selection(x)
  UpperProblemTime = @elapsed x.upper_problem!(x, CurrentNode)
  if x.current_upper_info.feasibility
    x.feasible_solution_found = true
    x.first_solution_node = x.maximum_node_id
    x.solution_value = x.current_upper_info.value
    x.continuous_solution[:] = x.current_upper_info.solution
    x.termination_status_code = MOI.LOCALLY_SOLVED
    x.result_status_code = MOI.FEASIBLE_POINT
  else
    x.feasible_solution_found = false
    x.first_solution_node = x.maximum_node_id
    x.termination_status_code = MOI.LOCALLY_INFEASIBLE
    x.result_status_code = MOI.INFEASIBLE_POINT
  end
    x.history.upper_time[1] = UpperProblemTime
end
