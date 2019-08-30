"""
    continuous_relative_bisect

Bisect node `S` on the maximum relative width dimension and return the two
resulting nodes.
"""
function continuous_relative_bisect(B::Optimizer, S::NodeBB)
  var_number = B.variable_number
  lower_variable_bounds = S.lower_variable_bounds
  upper_variable_bounds = S.upper_variable_bounds
  Pos = 0; Max = -Inf; TempMax = 0.0
  for i in 1:var_number
    if (~B.fixed_variable[i]) && (B.bisection_variable[i])
      TempMax = (upper_variable_bounds[i] - lower_variable_bounds[i])/(B.variable_info[i].upper_bound - B.variable_info[i].lower_bound)
      if TempMax > Max
        Pos = i; Max = TempMax
      end
    end
  end
  branch_pnt::Float64 = B.mid_cvx_factor*B.current_lower_info.solution[Pos] +
                        (1.0-B.mid_cvx_factor)*(lower_variable_bounds[Pos]+upper_variable_bounds[Pos])/2.0
  N1::Interval{Float64} = Interval{Float64}(lower_variable_bounds[Pos], branch_pnt)
  N2::Interval{Float64} = Interval{Float64}(branch_pnt, upper_variable_bounds[Pos])
  S.lower_bound = max(S.lower_bound, B.current_lower_info.value)
  S.upper_bound = min(S.upper_bound, B.current_upper_info.value)
  lvb_1 = Float64[ (i == Pos) ? N1.lo : lower_variable_bounds[i] for i in 1:var_number]
  uvb_1 = Float64[ (i == Pos) ? N1.hi : upper_variable_bounds[i] for i in 1:var_number]
  lvb_2 = Float64[ (i == Pos) ? N2.lo : lower_variable_bounds[i] for i in 1:var_number]
  uvb_2 = Float64[ (i == Pos) ? N2.hi : upper_variable_bounds[i] for i in 1:var_number]
  X1 = NodeBB(lvb_1, uvb_1, S.lower_bound, S.upper_bound, S.depth + 1, -1, false)
  X2 = NodeBB(lvb_2, uvb_2, S.lower_bound, S.upper_bound, S.depth + 1, -1, false)
  return X1, X2
end

#=
TO DO: ADD CODE FOR PSEUDO-COST BRANCHING
function pseudo_cost_branch(B::Optimizer,S::NodeBB)
    indx, psval = -1,-Inf
    for var in B.bisection_variable
      LowerCount = B.ProbCountLower[var]
      CostLower =(LowerCount > 0.0) ? (B.PseudoCostLower[var]/LowerCount) : 0.0
      UpperCount = B.ProbCountLower[var]
      CostUpper =(UpperCount > 0.0) ? (B.PseudoCostUpper[var]/UpperCount) : 0.0
      Score = max(CostLower,B.PseudoTol)*min(CostUpper,B.PseudoTol)
      if (Score > psval)
          psval = Score
          indx = var
      end
    end
    B.ProbCountLower[indx] += 1; B.ProbCountLower[var]
    return index, psval
end

# hybrid reliability rule branching on integer variables,
function PseudoCostBisect(B::Optimizer,S::NodeBB)
  # if fractional integer terms --> branch on them (hybrid reliability rule)
  FractionalInteger = is_integer_feasible(B)
  Pos,Value = PseudoCostBranch(B,S)
  if FractionalInteger
      S.LowerBound = max(S.LowerBound, B.CurrentLowerInfo.Value)
      S.UpperBound = min(S.UpperBound, B.CurrentUpperInfo.Value)
      X1 = copy(S); X1.LowerVars[Pos] = 0.0; X1.UpperVars[Pos] = 0.0; X1.Depth += 1; X1.LastBranch = Pos; X1.DirBranch = false
      X2 = copy(S); X2.LowerVars[Pos] = 1.0; X2.UpperVars[Pos] = 1.0; X2.Depth += 1; X2.LastBranch = Pos; X2.DirBranch = true
  else
      N1::Interval{Float64},N2::Interval{Float64} = bisect(Interval(S.LowerVars[Pos],S.UpperVars[Pos]),Value)
      S.LowerBound = max(S.LowerBound, B.CurrentLowerInfo.Value)
      S.UpperBound = min(S.UpperBound, B.CurrentUpperInfo.Value)
      X1 = copy(S); X1.LowerVars[Pos] = Nl.lo; X1.UpperVars[Pos] = Nl.hi; X1.Depth += 1; X1.LastBranch = Pos; X1.DirBranch = false
      X2 = copy(S); X2.LowerVars[Pos] = N2.lo; X2.UpperVars[Pos] = N2.hi; X2.Depth += 1; X2.LastBranch = Pos; X2.DirBranch = true
  end
end
=#
