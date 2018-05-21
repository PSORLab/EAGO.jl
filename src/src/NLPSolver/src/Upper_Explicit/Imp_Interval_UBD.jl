function Imp_Interval_UBD(X::Vector{V},k::Int,pos::Int,opt,temp) where {V<:AbstractInterval}
      nx::Int64 = opt[1].Imp_nx
      np::Int64 = opt[1].numVar - nx

      Eflag::Bool = false
      Iflag::Bool = false
      eDflag::Bool = false
      pnt::Vector{Float64} = mid.(X[(opt[1].Imp_nx+1):(opt[1].numVar)])
      PmidIntv = [V(pnt[i],pnt[i]) for i=1:np]
      Xc,Xc1,Eflg,Iflg,eDflg,incLow,incHigh = Param_Intv_Contractor(opt[1].Imp_h,opt[1].Imp_hj,
                                                                   X[1:opt[1].Imp_nx],PmidIntv,
                                                                   Eflag,Iflag,eDflag,
                                                                   opt[1].solver.PIntOpt)
      # solve optimization problem via interval extension
      feas = true
      FInt::V = opt[1].Imp_f(Xc,PmidIntv)
      val::Float64 = FInt.hi
      (opt[1].Imp_nCons > 0) && error("Problem must be unconstrained (aside from state variable)for interval-midpoint upper bound.")
      return val, vcat(mid.(Xc),pnt), feas, Any[feas,val]
end
