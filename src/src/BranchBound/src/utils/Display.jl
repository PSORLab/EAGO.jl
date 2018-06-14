"""
    EAGO.print_sol!(x::BnBSolver,y::BnBModel,ubdcnt::Int64,lbdcnt::Int64,ubdtime::Float64,lbdtime::Float64)

Prints solution information for the B&B problem. Displays first node found, solution value,
solution, and time spent solving subproblems.
"""
function print_sol!(x::BnBSolver,y::BnBModel,
                           ubdcnt::Int64,lbdcnt::Int64,
                           ubdtime::Float64,lbdtime::Float64)
  temp = 0.0
  println("First Solution Found at Node $(y.first_num)")
  if (x.Verbosity=="Normal"||x.Verbosity=="Full")
    println("UBD = $(y.soln_val)")
    println("Solution is :")
    if (y.feas_fnd)
      for i=1:length(y.Init_Box)
        temp = y.soln[i]
        println("    X[$i] = $temp")
      end
    end
    println("Total LBD problems solved = $lbdcnt in $lbdtime seconds.")
    println("Total UBD problems solved = $ubdcnt in $ubdtime seconds.")
    println("Total time spent preprocessing =  $(y.Pretime[end]) seconds.")
    println("Total time spent postprocessing = $(y.Posttime[end]) seconds.")
  end
end

"""
    EAGO.print_node!(x::BnBSolver,id::Int64,lbd::Float64,box::Vector{Interval{V}}) where {V<:AbstractFloat}

Prints node information for the B&B problem. Node id, bound, and interval box.
"""
function print_node!(x::BnBSolver,id::Int64,lbd::Float64,box::Vector{Interval{V}}) where {V<:AbstractFloat}
  if (x.Verbosity == "Full")
    println("Node ID: $(id), Lower Bound: $(lbd), IntervalBox: $(box)")
  end
end

function print_node!(x::BnBSolver,id::Int64,lbd::Float64,
                            box::Vector{MCInterval{V}}) where {V<:AbstractFloat}
  if (x.Verbosity == "Full")
    println("Node ID: $(id), Lower Bound: $(lbd), IntervalBox: $(box)")
  end
end

"""
    EAGO.print_int!(B::BnBSolver,k_int::Int64,k_nod::Int64,nid::Int64,lbdp::Float64,lbd::Float64,ubd::Float64,feasL::Bool,feasU::Bool)

Prints the iteration information if the Verbosity is set to "Normal" or "Full".
The header is displayed every hdr_intv, the iteration info is displayed every
itr_intv
"""
function print_int!(B::BnBSolver,k_int::Int64,k_nod::Int64,
                           nid::Int64,lbdp::Float64,lbd::Float64,ubd::Float64,feasL::Bool,feasU::Bool)
  if ((B.Verbosity == "Full")||(B.Verbosity == "Normal"))
    # prints header line every B.hdr_intv times
    if (mod(k_int,B.hdr_intv)==Int64(0)||k_int==Int64(1))
      println("Iteration   NodeID    Current_LBD     Global_LBD     Global_UBD      NodesLeft     Absolute_Gap    Absolute_Ratio     LBD_Feas     UBD_Feas")
    end
    # prints iteration summary every B.itr_intv times
    sbool1 = feasL ? "true" : "false"
    sbool2 = feasU ? "true" : "false"
    if ((mod(k_int,B.itr_intv)==Int64(0)))
      ptr_arr_temp = [k_int nid lbdp lbd ubd k_nod abs(ubd-lbd) abs(ubd-lbd)/(min(abs(lbd),abs(ubd))) sbool1 sbool2]
      ptr_arr1 = join([@sprintf("%6u",x) for x in ptr_arr_temp[1:2]], ",   ")
      ptr_arr2 = join([@sprintf("%3.7f",x) for x in ptr_arr_temp[3:5]], ",     ")
      ptr_arr3 = join([@sprintf("%6u",x) for x in ptr_arr_temp[6:6]], ",")
      ptr_arr4 = join([@sprintf("%3.7f",x) for x in ptr_arr_temp[7:8]], ",       ")
      ptr_arr5 = join([@sprintf("%s",x) for x in ptr_arr_temp[9:10]], ",       ")
      println(string(ptr_arr1,",      ",ptr_arr2,",      ",ptr_arr3),",        ", ptr_arr4,",        ", ptr_arr5)
    end
  end
end

"""
    EAGO.print_results!(B::BnBSolver,sol::Float64,pnt,feas::Bool,lbd_bool::Bool)

Prints the results of a single bounding problem.
"""
function print_results!(B::BnBSolver,sol::Float64,pnt,feas::Bool,lbd_bool::Bool)
  if (B.Verbosity == "Full")
    if (lbd_bool)
      println("Lower Bound: $(sol), Solution: $(pnt), Feasibility: $(feas)")
    else
      println("Upper Bound: $(sol), Solution: $(pnt), Feasibility: $(feas)")
    end
  end
end
