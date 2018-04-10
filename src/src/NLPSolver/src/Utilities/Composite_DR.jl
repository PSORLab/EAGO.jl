"""
    repeat_DR(X::Vector{Interval{Float64}},X0::Vector{Interval{Float64}},
              opt,k::Int64,rep::Bool)

Checks to see how much volume of the box was reduced and if it was reduced below
tolerance `tol_reduce_rept`.
"""
function repeat_DR(X::Vector{Interval{Float64}},
                   X0::Vector{Interval{Float64}},
                   opt,
                   k::Int64,
                   rep::Bool)

         if (opt.solver.max_reduce_rept < k)
           vol_diff::Float64 = 0.0
           for i=1:opt.numVar
             vol_diff += diam(X[i])/diam(X0[i])
           end
           if (opt.solver.tol_reduce_rept>vol_diff)
             rep = true
           end
         end
         return rep
end

function repeat_DR(X::Vector{MCInterval{Float64}},
                   X0::Vector{MCInterval{Float64}},
                   opt,
                   k::Int64,
                   rep::Bool)

         if (opt.solver.max_reduce_rept < k)
           vol_diff::Float64 = 0.0
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
                          bnbs::BnBSolver,bnbm::BnBModel{T}) where {T}
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
  println("Finished DAG Contractor")
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
  println("Finished Standard Range Reduction")
  # Implicit Range Reduction
  if (opt[1].solver.Implicit_Options.flag)
    if (feas == true && (opt[1].solver.Implicit_Options.Imp_RR_depth >= pos))
      Imp_Linear_RR!(X,opt,UBD)
      for i=(opt[1].solver.Implicit_Options.nx+1):(opt[1].numVar)
        if (isempty(X[i]))
          feas = false
          break
        end
      end
    end
  end
  println("Finished Implicit Range Reduction")
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
  println("Finished Standard Probing")
  # Implicit Probing
  if (opt[1].solver.Implicit_Options.flag)
    if (feas == true && (opt[1].solver.Implicit_Options.Imp_probe_depth >= pos))
      #Imp_LP_Probe!(X,opt,UBD)
      for i=(opt[1].solver.Implicit_Options.nx+1):(opt[1].numVar)
        if (isempty(X[i]))
          feas = false
          break
        end
      end
    end
  end
  println("Finished Implicit Probing")
  # Implicit Interval Contractor
  if (opt[1].solver.Implicit_Options.flag && (feas == true))
    Eflag = false
    Iflag = false
    eDflag = false
    if (opt[1].solver.Implicit_Options.Inplace)
      if (opt[1].solver.Implicit_Options.Intv_Cntr == "NewtonGS")
        #println("Started Inplace Newton")
        newtonGS = PId_NewtonGS(X[1:opt[1].solver.Implicit_Options.nx],
                                X[(opt[1].solver.Implicit_Options.nx+1):(opt[1].numVar)],
                                opt[1].solver.Implicit_Options.hj,
                                opt[1].solver.Implicit_Options.h,
                                opt[1].solver.Implicit_Options.ParamInt,
                                Eflag,Iflag,eDflag)
        X[1:opt[1].solver.Implicit_Options.nx] = newtonGS[1]
        if (newtonGS[3] == true)
          feas = false
        else (newtonGS[5] == true)
          feas = false
          # Solves LBD & UBD for node #1
          LBD1, LBDsol1, LBDfeas1, tempL1 = bnbs.Lower_Prob(X,k,pos,bnbs.opt,bnbm.UBDg)
          if (LBDfeas1)
            UBD1, UBDsol1, UBDfeas1, tempU1 = bnbs.Upper_Prob(X,k,pos,bnbs.opt,bnbm.UBDg)
          else
            UBD1 = copy(UBDn)
            UBDsol1 = mid.(X)
            UBDfeas1 = false
            tempU1 = []
          end

          # Solves LBD & UBD for node #1
          LBD2, LBDsol2, LBDfeas2, tempL2 = bnbs.Lower_Prob(X,k,pos,bnbs.opt,bnbm.UBDg)
          if (LBDfeas1)
            UBD2, UBDsol2, UBDfeas2, tempU2 = bnbs.Upper_Prob(X,k,pos,bnbs.opt,bnbm.UBDg)
          else
            UBD2 = copy(UBDn)
            UBDsol2 = mid.(X)
            UBDfeas2 = false
            tempU2 = []
          end

          # Stores if feasible
          if ((LBDfeas1 && UBDfeas1) && (LBDfeas2 && UBDfeas2))
            push!(bnbm.box,vcat(newtonGS[1],X[(opt[1].solver.Implicit_Options.nx+1):(opt[1].numVar)])
                          ,vcat(newtonGS[2],X[(opt[1].solver.Implicit_Options.nx+1):(opt[1].numVar)]))
            push!(bnbm.LBD,LBD1,LBD2)
            push!(bnbm.UBD,UBD1,UBD2)
            push!(bnbm.id,bnbm.max_id+1,bnbm.max_id+2)
            push!(bnbm.pos,pos+1,pos+1)
            bnbm.max_id += 2
          elseif (LBDfeas1 && UBDfeas1)
            push!(bnbm.box,vcat(newtonGS[1],X[(opt[1].solver.Implicit_Options.nx+1):(opt[1].numVar)]))
            push!(bnbm.LBD,LBD1)
            push!(bnbm.UBD,UBD1)
            push!(bnbm.id,bnbm.max_id+1)
            push!(bnbm.pos,pos+1)
            bnbm.max_id += 1
          elseif (LBDfeas2 && UBDfeas2)
            push!(bnbm.box,vcat(newtonGS[2],X[(opt[1].solver.Implicit_Options.nx+1):(opt[1].numVar)]))
            push!(bnbm.LBD,LBD2)
            push!(bnbm.UBD,UBD2)
            push!(bnbm.id,bnbm.max_id+1)
            push!(bnbm.pos,pos+1)
            bnbm.max_id += 1
          end
        end
      elseif (opt[1].solver.Implicit_Options.Intv_Cntr == "KrawczykCW")
        println("Started Inplace Krawczyk")
        krawczykCW = PId_KrawczykCW(X[1:opt[1].solver.Implicit_Options.nx],
                              X[(opt[1].solver.Implicit_Options.nx+1):(opt[1].numVar)],
                              opt[1].solver.Implicit_Options.hj,
                              opt[1].solver.Implicit_Options.h,
                              opt[1].solver.Implicit_Options.ParamInt,
                              Eflag,Iflag)
        X[1:opt[1].solver.Implicit_Options.nx] = krawczykCW[1]
        if (krawczykCW[2] == true)
          feas = false
        end
      end
    else
      if (opt[1].solver.Implicit_Options.Intv_Cntr == "NewtonGS")
        println("Started Outplace Newton")
        newtonGS = PI_NewtonGS(X[1:opt[1].solver.Implicit_Options.nx],
                                X[(opt[1].solver.Implicit_Options.nx+1):(opt[1].numVar)],
                                opt[1].solver.Implicit_Options.hj,
                                opt[1].solver.Implicit_Options.h,
                                opt[1].solver.Implicit_Options.ParamInt,
                                Eflag,Iflag,eDflag)
        X[1:opt[1].solver.Implicit_Options.nx] = newtonGS[1]
        if (newtonGS[3] == true)
          feas = false
        else (newtonGS[5] == true)
          feas = false
          # Solves LBD & UBD for node #1
          LBD1, LBDsol1, LBDfeas1, tempL1 = bnbs.Lower_Prob(X,k,pos,bnbs.opt,bnbm.UBDg)
          if (LBDfeas1)
            UBD1, UBDsol1, UBDfeas1, tempU1 = bnbs.Upper_Prob(X,k,pos,bnbs.opt,bnbm.UBDg)
          else
            UBD1 = copy(UBDn)
            UBDsol1 = mid.(X)
            UBDfeas1 = false
            tempU1 = []
          end

          # Solves LBD & UBD for node #1
          LBD2, LBDsol2, LBDfeas2, tempL2 = bnbs.Lower_Prob(X,k,pos,bnbs.opt,bnbm.UBDg)
          if (LBDfeas1)
            UBD2, UBDsol2, UBDfeas2, tempU2 = bnbs.Upper_Prob(X,k,pos,bnbs.opt,bnbm.UBDg)
          else
            UBD2 = copy(UBDn)
            UBDsol2 = mid.(X)
            UBDfeas2 = false
            tempU2 = []
          end

          # Stores if feasible
          if ((LBDfeas1 && UBDfeas1) && (LBDfeas2 && UBDfeas2))
            push!(bnbm.box,vcat(newtonGS[1],X[(opt[1].solver.Implicit_Options.nx+1):(opt[1].numVar)])
                          ,vcat(newtonGS[2],X[(opt[1].solver.Implicit_Options.nx+1):(opt[1].numVar)]))
            push!(bnbm.LBD,LBD1,LBD2)
            push!(bnbm.UBD,UBD1,UBD2)
            push!(bnbm.id,bnbm.max_id+1,bnbm.max_id+2)
            push!(bnbm.pos,pos+1,pos+1)
            bnbm.max_id += 2
          elseif (LBDfeas1 && UBDfeas1)
            push!(bnbm.box,vcat(newtonGS[1],X[(opt[1].solver.Implicit_Options.nx+1):(opt[1].numVar)]))
            push!(bnbm.LBD,LBD1)
            push!(bnbm.UBD,UBD1)
            push!(bnbm.id,bnbm.max_id+1)
            push!(bnbm.pos,pos+1)
            bnbm.max_id += 1
          elseif (LBDfeas2 && UBDfeas2)
            push!(bnbm.box,vcat(newtonGS[2],X[(opt[1].solver.Implicit_Options.nx+1):(opt[1].numVar)]))
            push!(bnbm.LBD,LBD2)
            push!(bnbm.UBD,UBD2)
            push!(bnbm.id,bnbm.max_id+1)
            push!(bnbm.pos,pos+1)
            bnbm.max_id += 1
          end
        end
      elseif (opt[1].solver.Implicit_Options.Intv_Cntr == "KrawczykCW")
        println("Started Outplace Krawczyk")
        krawczykCW = PI_KrawczykCW(X[1:opt[1].solver.Implicit_Options.nx],
                              X[(opt[1].solver.Implicit_Options.nx+1):(opt[1].numVar)],
                              opt[1].solver.Implicit_Options.hj,
                              opt[1].solver.Implicit_Options.h,
                              opt[1].solver.Implicit_Options.ParamInt,
                              Eflag,Iflag)
        X[1:opt[1].solver.Implicit_Options.nx] = krawczykCW[1]
        if (krawczykCW[2] == true)
          feas = false
        end
      end
    end
  end
  println("Finished Interval Contractor")
  #println("feas: $feas")
  #println("end preprocessing: ")

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
function composite_DR_post(feas_Post::Bool,X::Vector{Interval{Float64}},k::Int64,pos::Int64,
                      opt,temp_objL,temp_objU,LBDg,UBDg)
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

      if (opt[1].solver.Implicit_Options.flag)
        if (opt[1].solver.variable_depth>=pos)
          Variable_DR_Imp!(X,dual_lo,dual_hi,LBD,UBD,opt[1].solver.Implicit_Options.nx)
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

function composite_DR_post(feas_Post::Bool,X::Vector{MCInterval{Float64}},k::Int64,pos::Int64,
                      opt,temp_objL,temp_objU,LBDg,UBDg)
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

      if (opt[1].solver.Implicit_Options.flag)
        if (opt[1].solver.variable_depth>=pos)
          Variable_DR_Imp!(X,dual_lo,dual_hi,LBD,UBD,opt[1].solver.Implicit_Options.nx)
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
