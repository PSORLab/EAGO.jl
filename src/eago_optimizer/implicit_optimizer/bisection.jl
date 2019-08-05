"""
    implicit_bisection

Bisect node `S` on the maximum relative width dimension of the first 1:np dimensions
and return the two resulting nodes.
"""
function implicit_bisection(B::Optimizer,S::NodeBB)
  Pos = 0; Max = -Inf; TempMax = 0.0
  nx = num_state_variables(B.working_evaluator_block.evaluator)
  np = num_decision_variables(B.working_evaluator_block.evaluator)
  for i in 1:np
    shift = nx + i
    if (~B.fixed_variable[shift]) && (B.bisection_variable[shift])
      TempMax = (S.upper_variable_bounds[shift] - S.lower_variable_bounds[shift])/(B.variable_info[shift].upper_bound - B.variable_info[shift].lower_bound)
      if TempMax > Max
        Pos = shift; Max = TempMax
      end
    end
  end
  CutInterval = IntervalType(S.lower_variable_bounds[Pos],S.upper_variable_bounds[Pos])
  N1::IntervalType, N2::IntervalType = bisect(CutInterval)
  S.lower_bound = max(S.lower_bound, B.current_lower_info.value)
  S.upper_bound = min(S.upper_bound, B.current_upper_info.value)
  X1 = NodeBB(S.lower_variable_bounds, S.upper_variable_bounds, S.lower_bound, S.upper_bound, S.depth + 1, -1, false)
  X2 = deepcopy(X1)
  X1.lower_variable_bounds[Pos] = N1.lo; X1.upper_variable_bounds[Pos] = N1.hi
  X2.lower_variable_bounds[Pos] = N2.lo; X2.upper_variable_bounds[Pos] = N2.hi
  return X1, X2
end
