"""
    Implicit_SIP_Solve

Solves a semi-infinite program subject to imbedded equality constraint that
defines an implicit function via the algorithm presented in Stuber2015 using
the EAGOGlobalSolver to solve the lower bounding problem, lower level problem,
and the upper bounding problem. The options for the algorithm and the global
solvers utilized are set by manipulating a SIPopt containing the options info.
Inputs:
* `f::Function`: Objective in the decision and state space variable.
                 Takes arguments vector, `(x,y)`, The function must be untyped.
* `h::Function`: Equality constraints on decision and uncertainty
                 variables. The arguements of `h` are `(out,x,y,p)` where `x` is the
                 control variable, `y` is the state variable, and p is the
                 uncertain varable.
* `hj::Function`: Jacobian of `h(out,x,y,p)` with respect to control & uncertain
                  variables.
* `gSIP::Function`: Semi-infinite constraint. The arguements of `g` are `(x,y,p)`
                    where `u` is the control variable, `x` is the state variable,
                    and `p` is the uncertain varable. The function must be untyped.
* `X::Vector{Interval{Float64}}`: Box constraints for state variables
* `Y::Vector{Interval{Float64}}`: Box constraints for decision variables
* `P::Vector{Interval{Float64}}`: Box constraints for uncertainty variables
* `SIPopt::SIP_opts`: Option type containing problem information

Returns:
A SIP_result composite type containing solution information.
"""
# Currently DAG contractor
function Implicit_SIP_Solve(f,h,hj,gSIP,X,Y,P,SIPopt::SIP_opts)

    # initializes solution
    UBDg = Inf
    LBDg = -Inf
    k = 0
    P_LBD = SIPopt.P_LBD
    P_UBD = SIPopt.P_UBD
    np = length(P)
    nx = length(X)
    ny = length(Y)
    P_low = [P[i].lo for i=1:np]
    P_high = [P[i].hi for i=1:np]
    Y_low = [Y[i].lo for i=1:nx]
    Y_high = [Y[i].hi for i=1:nx]
    X_low = [X[i].lo for i=1:nx]
    X_high = [X[i].hi for i=1:nx]
    pbar = mid.(P)
    xbar = mid.(X)
    INNg1 = Inf
    INNg2 = Inf
    feas = true

    # Turns implicit solver routines on
    SIPopt.LBP_Opt.ImplicitFlag = true
    SIPopt.LLP_Opt.ImplicitFlag = true
    SIPopt.UBP_Opt.ImplicitFlag = true

    # Sets number of state variables
    SIPopt.LBP_Opt.PSmcOpt.nx = nx
    SIPopt.LBP_Opt.PIntOpt.nx = nx
    SIPopt.LLP_Opt.PSmcOpt.nx = nx
    SIPopt.LLP_Opt.PIntOpt.nx = nx
    SIPopt.UBP_Opt.PSmcOpt.nx = nx
    SIPopt.UBP_Opt.PIntOpt.nx = nx

    # creates results file to store output
    sip_sto = SIP_result()

    # checks inputs
    if (SIPopt.r0<=1)
      error("r0 must be greater than 1")
    elseif (SIPopt.eps_g0<=0)
      error("eps_g must be greater than 0")
    else
      eps_g = SIPopt.eps_g0
      r = SIPopt.r0
    end

    ##### checks for convergence #####
    for k=1:SIPopt.kmax

      ##### check for termination #####
      if (abs(UBDg-LBDg)<SIPopt.tol)
        println("Algorithm Converged")
        break
      end
      println("ran to here 1")
      ##### lower bounding problem #####
      LBD_X_low,LBD_X_high,refnx,snx = Reform_Imp_Y(X,Y,P_LBD)
      gL_LBP = [-Inf for i=1:((1+2*ny)*length(P_LBD))]
      gU_LBP = [0.0 for i=1:((1+2*ny)*length(P_LBD))]
      mLBP = deepcopy(MathProgBase.NonlinearModel(SIPopt.LBP_Opt))
      mLBP.Opts.Imp_nx = refnx
      gLBP = x -> Reform_Imp_HG(gSIP,h,x[(refnx+1):end],x[1:refnx],pLBP,ny,0.0)
      println("snx: $snx")
      println("refnx: $refnx")
      println("LBD_X_low: $(length(LBD_X_low))")
      println("LBD_X_high: $(length(LBD_X_high))")
      println("gL_LBP: $(length(gL_LBP))")
      println("gU_LBP: $(length(gU_LBP))")
      MathProgBase.loadproblem!(mLBP, refnx, snx, LBD_X_low, LBD_X_high,
                                gL_LBP, gU_LBP, :Min, f, gLBP)
      mLBP.Opts.Imp_f = (y,p) -> f(p,y)
      mLBP.Opts.Imp_g = (y,p) -> Reform_Imp_G(g,p,y,P_UBD,ny,0.0)
      mLBP.Opts.Imp_h = (y,p) -> Reform_Imp_H(p,y,P_LBD,ny)
      mLBP.Opts.Imp_hj = (y,p) -> Reform_Imp_HJ(p,y,P_LBD,ny)

      if SIPopt.LBP_Opt.DAG_depth>0
        if (SIPopt.gSIPExp == Expr[])
          error("Must provide expression for gSIP in order to use DAG contractor")
        else
          mLBP.Opts.DAG_tlist = Generate_Fixed_TapeList(x->gSIP(x[1:nx],x[(nx+1):(nx+np)]),nx,gL_LBP,gU_LBP,P_LBD)
        end
      end

      MathProgBase.optimize!(mLBP)
      LBDg = MathProgBase.getobjval(mLBP)
      xbar = MathProgBase.getsolution(mLBP)
      feas = getfeasibility(mLBP)
      tLBP = MathProgBase.getsolvetime(mLBP)
      sip_sto.LBP_time += tLBP
      sip_sto.LBD = LBDg
      sip_sto.xbar = xbar
      if (SIPopt.Verbosity == "Full" || SIPopt.Verbosity == "Normal")
        println("solved LBD: ",LBDg," ",xbar," ",feas)
      end
      if (~feas)
        println("Lower Bounding Problem Not Feasible. Algorithm Terminated")
        sip_sto.feas = false
        return sip_sto
      end

      ##### inner program #####
      mLLP1 = deepcopy(MathProgBase.NonlinearModel(SIPopt.LLP_Opt))
      MathProgBase.loadproblem!(mLLP1, np, 0, P_low, P_high,
                                Float64[], Float64[], :Min, p -> -gSIP(xbar,p), [])
      MathProgBase.optimize!(mLLP1)
      INNg1 = MathProgBase.getobjval(mLLP1)
      pbar = MathProgBase.getsolution(mLLP1)
      feas = getfeasibility(mLLP1)
      tLLP = MathProgBase.getsolvetime(mLLP1)
      INNg1 = - INNg1
      sip_sto.LLP_time += tLLP
      if (SIPopt.Verbosity == "Full" || SIPopt.Verbosity == "Normal")
        println("solved INN #1: ",INNg1," ",pbar," ",feas)
      end
      if (INNg1+SIPopt.inn_tol<=0)
        xstar = xbar
        UBDg = LBDg
        return LBDg,UBDg
      else
        push!(P_LBD,pbar)
      end
      ##### upper bounding problem #####
      UBD_X_low,UBD_X_nx,refnx,snx = Reform_Imp_Y(X,Y,P_UBD)
      gL_UBP = [-Inf for i=1:length(P_UBD)]
      gU_UBP = [0.0 for i=1:length(P_UBD)]
      gUBP = x -> Reform_Imp_G(gSIP,x[(refnx+1):end],x[1:refnx],pUBD,ny,eps_g)
      mUBP = deepcopy(MathProgBase.NonlinearModel(SIPopt.UBP_Opt))
      mUBP.Opts.Imp_nx = refnx
      mUBP.Opts.Imp_f = (y,p) -> f(p,y)
      mUBP.Opts.Imp_g = (y,p) -> Reform_Imp_G(gSIP,p,y,P_UBD,ny,eps_g)
      mUBP.Opts.Imp_h = (y,p) -> Reform_Imp_H(p,y,P_UBD,ny)
      mUBP.Opts.Imp_hj = (y,p) -> Reform_Imp_HJ(p,y,P_UBD,ny)
      MathProgBase.loadproblem!(mUBP, snx, length(P_UBD), UBD_X_low, UBD_X_nx,
                                gL_UBP, gU_UBP, :Min, f, gUBP)
      if SIPopt.UBP_Opt.DAG_depth>0
        if (SIPopt.gSIPExp == Expr[])
          error("Must provide expression for gSIP in order to use DAG contractor")
        else
          mUBP.Opts.DAG_tlist = Generate_Fixed_TapeList(x->gSIP(x[1:nx],x[(nx+1):(nx+np)]),nx,gL_UBP,gU_UBP,P_UBD)
        end
      end
      MathProgBase.optimize!(mUBP)
      UBD_temp = MathProgBase.getobjval(mUBP)
      xbar = MathProgBase.getsolution(mUBP)
      feas = getfeasibility(mUBP)
      tUBP = MathProgBase.getsolvetime(mUBP)
      sip_sto.UBP_time += tUBP
      sip_sto.UBD = UBD_temp
      sip_sto.xbar = xbar

      if (SIPopt.Verbosity == "Full" || SIPopt.Verbosity == "Normal")
        println("solved UBD: ",UBD_temp," ",xbar," ",feas)
      end
      if (feas)
        ##### inner program #####
        mLLP2 = deepcopy(MathProgBase.NonlinearModel(SIPopt.LLP_Opt))
        MathProgBase.loadproblem!(mLLP2, np, 0, P_low, P_high,
                                  Float64[], Float64[], :Min, p -> -gSIP(xbar,p), [])
        MathProgBase.optimize!(mLLP2)
        INNg2 = MathProgBase.getobjval(mLLP2)
        pbar = MathProgBase.getsolution(mLLP2)
        feas = getfeasibility(mLLP2)
        tLLP = MathProgBase.getsolvetime(mLLP2)
        sip_sto.LLP_time += tLLP
        INNg2 = - INNg2
        if (SIPopt.Verbosity == "Full" || SIPopt.Verbosity == "Normal")
          println("solved INN #2: ",INNg2," ",pbar," ",feas)
        end
        if (INNg2+SIPopt.inn_tol<0)
          if (UBD_temp <= UBDg)
            UBDg = UBD_temp
            xstar = xbar
          end
          eps_g = eps_g/r
        else
          push!(P_UBD,pbar)
        end
      else
        eps_g = eps_g/r
      end

      print_int!(SIPopt,k,LBDg,UBDg,eps_g,r)
      sip_sto.k += k
    end

    return sip_sto
  end
