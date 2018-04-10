"""
    Term_Check(x::BnBSolver,y::BnBModel,k_int::Int64)

Checks if algorithm continues for the branch and bound problem. Checks that:
* Nodes are remaining on the stack.
* The iteration limit (`x.max_iter`) is not exceeded if iteration limit
flag (`x.iter_lim`) is true.
* Maximum nodes (`x.max_nodes`) on stack are not exceeded.
* Absolute global tolerance (`x.BnB_atol`) is not reached.
* Relative global tolerance (`x.BnB_rtol`) is not reached.
* Target upper bound (`x.target_upper`) is not reached.
Outputs a description of termination check if it returns false.
"""
function Term_Check(x::BnBSolver,y::BnBModel,k_int::Int64)
  t1 = length(y.LBD)>0
  t2 = (x.iter_lim ? k_int<x.max_iter : true)
  t3 = (length(y.LBD)<x.max_nodes)
  t4 = (y.UBDg-y.LBDg) > x.BnB_atol
  t5 = (y.UBDg-y.LBDg) > abs(y.LBDg)*x.BnB_rtol || ~(y.LBDg > -Inf)
  t6 = x.target_upper <= y.UBDg
  if t1 & t2 & t3 & t4 & t5 & t6
    return true
  else
    if (x.Verbosity=="Normal"||x.Verbosity=="Full")
      if ~(length(y.LBD)>0)
        if (y.first_num>0)
          println("Empty Stack: Exhaustive Search Finished")
        else
          println("Empty Stack: Infeasible")
        end
      elseif ~(length(y.LBD)<x.max_nodes)
        println("Node Limit Exceeded")
      elseif ~(x.iter_lim ? k_int<x.max_iter : true)
        println("Maximum Iteration Exceeded")
      else
        println("Convergence Tolerance Reached")
      end
    end
    return false
  end
end

"""
    Conv_Check(x::BnBSolver,ubd::Float64,lbd::Float64)

Checks if convergence tolerance is reach for the branch and bound problem.
Checks that:
* Absolute global tolerance (`x.BnB_atol`) is reached.
* Relative global tolerance (`x.BnB_rtol`) is reached.
Outputs a description of termination check if it returns false.
"""
function Conv_Check(x::BnBSolver,ubd::Float64,lbd::Float64)
  return ((abs(ubd-lbd) <= x.BnB_atol) || (abs(ubd-lbd)/(min(abs(lbd),abs(ubd))) <= x.BnB_rtol))
end

"""
    Repeat_Node_Default(x::BnBSolver,y::BnBModel,Xin::Vector{Interval{Float64}},
                                 Xout::Vector{Interval{Float64}})

Default check for repeating a node. Always returns false.
"""
function Repeat_Node_Default(x::BnBSolver,y::BnBModel{Interval{T}},
                             Xin::Vector{Interval{T}},
                             Xout::Vector{Interval{T}}) where {T<:AbstractFloat}
  return false
end

function Repeat_Node_Default(x::BnBSolver,y::BnBModel{MCInterval{T}},
                             Xin::Vector{MCInterval{T}},
                             Xout::Vector{MCInterval{T}}) where {T<:AbstractFloat}
  return false
end
