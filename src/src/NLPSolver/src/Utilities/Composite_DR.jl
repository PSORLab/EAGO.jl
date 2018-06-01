"""
    repeat_DR(X::Vector{Interval{Float64}},X0::Vector{Interval{Float64}},
              opt,k::Int64,rep::Bool)

Checks to see how much volume of the box was reduced and if it was reduced below
tolerance `tol_reduce_rept`.
"""
function repeat_DR(X::Vector{V},
                   X0::Vector{V},
                   opt,
                   k::Q,
                   rep::Bool) where {Q<:Integer,V<:AbstractInterval}

         if (opt.solver.max_reduce_rept < k)
           vol_diff = 0.0
           for i=1:opt.numVar
             vol_diff += diam(X[i])/diam(X0[i])
           end
           if (opt.solver.tol_reduce_rept>vol_diff)
             rep = true
           end
         end
         return rep
end


"""
    composite_DR_pre(feas::Bool,X::Vector{Interval{Float64}},UBD::Float64,
                     k::Int64,pos::Int64,opt)

A function that performs all pre-processing routines in EAGO if the appropriate
flag is selected in the solver options. The inputs are:
* `feas::Bool`: The input feasibility of problem.
* `X::Vector{Interval{Float64}}`: The input variable bounds.
* `UBD::Float64`: The current global bounds.
* `k::Int64`: The iteration number.
* `pos::Int64`: The position in the tree
* `opt`: The option file containing solver options.
Outputs the feasibility and the X.
"""
function composite_DR_pre(feas::Bool,X::Vector{T},UBD::Float64,
                          k::Int64,pos::Int64,opt,LBDn::Float64,UBDn::Float64,
                          bnbs::BnBSolver,bnbm::BnBModel{T}) where {T<:AbstractInterval}

  # Feasibility-Based Bound Tightening via DAG Constraint Propagation
  if (feas == true && (opt[1].solver.DAG_depth >= pos))
    DAGContractor!(X,opt[1].DAG_tlist,opt[1].solver.DAG_pass)
    for i=1:opt[1].numVar
      if (isempty(X[i]))
        feas = false
        break
      end
    end
  end
  (bnbs.Verbosity == "Full") && (println("Finished DAG Contractor"))
  # Standard Range Reduction
  if (feas == true && (opt[1].solver.STD_RR_depth >= pos))
    STD_Linear_RR!(X,opt,UBD)
    for i=1:opt[1].numVar
      if (isempty(X[i]))
        feas = false
        break
      end
    end
  end
  (bnbs.Verbosity == "Full") && (println("Finished Standard Range Reduction"))
  # Implicit Range Reduction
  if (opt[1].solver.ImplicitFlag)
    if (feas == true && (opt[1].solver.Imp_RR_depth >= pos))
      Imp_Linear_RR!(X,opt,UBD)
      for i=(opt[1].solver.Imp_nx+1):(opt[1].numVar)
        if (isempty(X[i]))
          feas = false
          break
        end
      end
    end
  end
  (bnbs.Verbosity == "Full") && (println("Finished Implicit Range Reduction"))
  # Standard  Probing
  if (feas == true && (opt[1].solver.probe_depth >= pos))
    STD_LP_Probe!(X,opt,UBD)
    for i=1:opt[1].numVar
      if (isempty(X[i]))
        feas = false
        break
      end
    end
  end
  (bnbs.Verbosity == "Full") && (println("Finished Standard Probing"))
  # Implicit Probing
  if (opt[1].solver.ImplicitFlag)
    if (feas == true && (opt[1].solver.Imp_probe_depth >= pos))
      #Imp_LP_Probe!(X,opt,UBD)
      for i=(opt[1].solver.Imp_nx+1):(opt[1].numVar)
        if (isempty(X[i]))
          feas = false
          break
        end
      end
    end
  end
  (bnbs.Verbosity == "Full") && (println("Finished Implicit Probing"))
  # Implicit Interval Contractor
  if (opt[1].solver.ImplicitFlag && (feas == true))
    #println("X[1:opt[1].Imp_nx]: $(X[1:opt[1].Imp_nx])")
    Eflag::Bool = false
    Iflag::Bool = false
    eDflag::Bool = false
    Y1,Y2,Eflag,Iflag,eDflag,incLow,incHigh = Param_Intv_Contractor(opt[1].Imp_h,opt[1].Imp_hj,
                                                                    X[1:opt[1].Imp_nx],
                                                                    X[(opt[1].Imp_nx+1):(opt[1].numVar)],
                                                                    Eflag,Iflag,eDflag,
                                                                    opt[1].solver.PIntOpt)
    #println("Y1: $Y1")
    #println("Eflag: $Eflag")
    #println("Iflag: $Iflag")
    #println("eDflag: $eDflag")
    if Eflag
      feas = false
    elseif eDflag
      feas = false

      # Solves LBD & UBD for node #1
      LBD1, LBDsol1, LBDfeas1, tempL1 = bnbs.Lower_Prob(Y1,k,pos,bnbs.opt,bnbm.UBDg)
      if (LBDfeas1)
        UBD1, UBDsol1, UBDfeas1, tempU1 = bnbs.Upper_Prob(Y1,k,pos,bnbs.opt,bnbm.UBDg)
      else
        UBD1 = copy(UBDn)
        UBDsol1 = mid.(Y)
        UBDfeas1 = false
        tempU1 = []
      end

      # Solves LBD & UBD for node #1
      LBD2, LBDsol2, LBDfeas2, tempL2 = bnbs.Lower_Prob(Y2,k,pos,bnbs.opt,bnbm.UBDg)
      if (LBDfeas1)
        UBD2, UBDsol2, UBDfeas2, tempU2 = bnbs.Upper_Prob(Y2,k,pos,bnbs.opt,bnbm.UBDg)
      else
        UBD2 = copy(UBDn)
        UBDsol2 = mid.(Y2)
        UBDfeas2 = false
        tempU2 = []
      end

      # Stores if feasible
      if ((LBDfeas1 && UBDfeas1) && (LBDfeas2 && UBDfeas2))
        push!(bnbm.box,vcat(Y1,X[(opt[1].Imp_nx+1):(opt[1].numVar)])
                      ,vcat(Y2,X[(opt[1].Imp_nx+1):(opt[1].numVar)]))
        push!(bnbm.LBD,LBD1,LBD2)
        push!(bnbm.UBD,UBD1,UBD2)
        push!(bnbm.id,bnbm.max_id+1,bnbm.max_id+2)
        push!(bnbm.pos,pos+1,pos+1)
        bnbm.max_id += 2
      elseif (LBDfeas1 && UBDfeas1)
        push!(bnbm.box,vcat(Y1,X[(opt[1].Imp_nx+1):(opt[1].numVar)]))
        push!(bnbm.LBD,LBD1)
        push!(bnbm.UBD,UBD1)
        push!(bnbm.id,bnbm.max_id+1)
        push!(bnbm.pos,pos+1)
        bnbm.max_id += 1
      elseif (LBDfeas2 && UBDfeas2)
        push!(bnbm.box,vcat(Y2,X[(opt[1].Imp_nx+1):(opt[1].numVar)]))
        push!(bnbm.LBD,LBD2)
        push!(bnbm.UBD,UBD2)
        push!(bnbm.id,bnbm.max_id+1)
        push!(bnbm.pos,pos+1)
        bnbm.max_id += 1
      end
    else
      X = vcat(Y1,X[(opt[1].Imp_nx+1):(opt[1].numVar)])
    end
  end
  (bnbs.Verbosity == "Full") && (println("Finished Implicit Contractor"))
  return feas,X
end

"""
    composite_DR_post(feas_Post::Bool,X::Vector{Interval{Float64}},k::Int64,pos::Int64,
                      opt,temp_objL,temp_objU,LBDg,UBDg)

A function that performs all post-processing routines in EAGO if the appropriate
flag is selected in the solver options. Supports both explicit and impliit solvers.
The inputs are:
* `feas_Post::Bool`: The input feasibility of problem.
* `X::Vector{Interval{Float64}}`: The input variable bounds.
* `k::Int64`: The iteration number.
* `pos::Int64`: The position in the tree
* `opt`: The option file containing solver options.
* `temp_objL`: The current global lower bounds.
* `temp_objU`: The current global upper bounds.
* `LBDg::Float64`: The current global lower bounds.
* `UBDg::Float64`: The current global upper bounds.
Outputs the feasibility and the X.
"""
function composite_DR_post(feas_Post::Bool,X::Vector{V},k::Q1,pos::Q2,
                      opt,temp_objL,temp_objU,LBDg,UBDg) where {Q1<:Integer,Q2<:Integer,V<:AbstractInterval}
  #println("start composite DR post")
  # Duality-Based Range Reduction
  if (feas_Post)
    if (length(temp_objL) > 0)
      dual_lo = temp_objL[1]
      dual_hi = temp_objL[2]
      LBD = temp_objL[3]
      if (temp_objU[1])
        UBD = min(UBDg,temp_objU[2])
      else
        UBD = UBDg
      end

      if (opt[1].solver.ImplicitFlag)
        if (opt[1].solver.variable_depth>=pos)
          Variable_DR_Imp!(X,dual_lo,dual_hi,LBD,UBD,opt[1].Imp_nx)
        end
      else
        if (opt[1].solver.variable_depth>=pos)
          Variable_DR!(X,dual_lo,dual_hi,LBD,UBD)
        end
      end
    end
  end
  #println("end composite DR post")
  return feas_Post, X
end
