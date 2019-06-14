"""
    default_repeat_check

Used to determine if optimizer should repeat processing due to node size reduction.
Currently returns false.
"""
function default_repeat_check(x::Optimizer,y::NodeBB)
  return false
end

"""
    default_termination_check

Checks for termination based on termination criteria.
"""
function default_termination_check(x::Optimizer)

    L = x.global_lower_bound
    U = x.global_upper_bound
    t1 = length(x.stack) > 0                                   # stack is non-empty
    t2 = x.current_iteration_count < x.iteration_limit            # maximum iterations not exceeded
    t3 = length(x.stack) < x.node_limit                         # maximum node number not exceeded
    t4 = (U - L) > x.absolute_tolerance                         # absolute tolerance satisfied
    t5 = (abs(U - L)/(min(abs(L),abs(U))) > x.relative_tolerance) || ~(L > -Inf)   # relative tolerance satisfied
    if t1 & t2 & t3 & t4 & t5
      return true
    else
      if ~t1
        if (x.first_solution_node > 0)
          x.termination_status_code = MOI.OPTIMAL
          x.result_status_code = MOI.FEASIBLE_POINT
          (x.verbosity >= 3) && println("Empty Stack: Exhaustive Search Finished")
        else
          x.termination_status_code = MOI.INFEASIBLE
          x.result_status_code = MOI.INFEASIBILITY_CERTIFICATE
          (x.verbosity >= 3) && println("Empty Stack: Infeasible")
        end
      elseif ~t3
        (x.verbosity >= 3) && println("Node Limit Exceeded")
        x.termination_status_code = MOI.NODE_LIMIT
        x.result_status_code = MOI.UNKNOWN_RESULT_STATUS
      elseif ~t2
        (x.verbosity >= 3) && println("Maximum Iteration Exceeded")
        x.termination_status_code = MOI.ITERATION_LIMIT
        x.result_status_code = MOI.UNKNOWN_RESULT_STATUS
      else
        x.termination_status_code = MOI.OPTIMAL
        x.result_status_code = MOI.FEASIBLE_POINT
        (x.verbosity >= 3) && println("Convergence Tolerance Reached")
      end
    end
    return false
end

"""
    default_convergence_check

Checks relative and absolute convergence.
"""
function default_convergence_check(x::Optimizer)
  L = x.current_lower_info.value
  U = x.current_upper_info.value
  t1 = (U - L) <= x.absolute_tolerance
  t2 = (abs(U - L)/(max(abs(L),abs(U))) <= x.relative_tolerance)
  return t1 || t2
end
