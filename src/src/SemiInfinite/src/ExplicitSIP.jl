"""
    Explicit_SIP_Solve

Solves a semi-infinite program via the algorithm presented in Mitsos2011 using
the EAGOGlobalSolver to solve the lower bounding problem, lower level problem,
and the upper bounding problem. The options for the algorithm and the global
solvers utilized are set by manipulating a SIPopt containing the options info.
Inputs:
* `f::Function`: Objective in the decision variable. Takes a single argument
                 vector that must be untyped.
* `gSIP::Function`: The semi-infinite constraint. Takes two arguments: the first
                    being a vector containing the decision variable and the
                    second being a vector containing the uncertainity
                    variables. The function must be untyped.
* `X::Vector{Interval}`: Box constraints for decision variables
* `P::Vector{Interval}`: Box constraints for uncertainty variables
* `SIPopt::SIP_opts`: Option type containing problem information
Returns: A SIP_result composite type containing solution information.
"""
function Explicit_SIP_Solve(f,gSIP,X,P,SIPopt::SIP_opts)

  # initializes solution
  UBDg = Inf
  LBDg = -Inf
  k = 0
  P_LBD = SIPopt.P_LBD
  P_UBD = SIPopt.P_UBD
  np = length(P)
  nx = length(X)
  P_low = [P[i].lo for i=1:np]
  P_high = [P[i].hi for i=1:np]
  X_low = [X[i].lo for i=1:nx]
  X_high = [X[i].hi for i=1:nx]
  pbar = mid.(P)
  xbar = mid.(X)
  INNg1 = Inf
  INNg2 = Inf
  feas = true

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

  #println("ran to algo start #1")
  ##### checks for convergence #####
  for k=1:SIPopt.kmax

    ##### check for termination #####
    if (abs(UBDg-LBDg)<SIPopt.tol)
      println("Algorithm Converged")
      break
    end

    ##### lower bounding problem #####
    gLBP = x -> BndProb_reform(x,nothing,gSIP,P_LBD,0.0) # reformulate constraints
    gL_LBP = [-Inf for i=1:length(P_LBD)]
    gU_LBP = [0.0 for i=1:length(P_LBD)]
    mLBP = deepcopy(MathProgBase.NonlinearModel(SIPopt.LBP_Opt))
    #  println("sarted lower problem load")
    MathProgBase.loadproblem!(mLBP, nx, length(P_LBD), X_low, X_high,
                              gL_LBP, gU_LBP, :Min, f, gLBP)

    #println("finished lower problem load")
    if SIPopt.LBP_Opt.DAG_depth>0
      if (SIPopt.gSIPExp == Expr[])
        error("Must provide expression for gSIP in order to use DAG contractor")
      else
        mLBP.Opts.DAG_tlist = Generate_Fixed_TapeList(x->gSIP(x[1:nx],x[(nx+1):(nx+np)]),nx,gL_LBP,gU_LBP,P_LBD)
      end
    end
    #println("lower tapelist")

    MathProgBase.optimize!(mLBP)
    #println("lower solve")
    LBDg = MathProgBase.getobjval(mLBP)
    xbar = MathProgBase.getsolution(mLBP)
    feas = getfeasibility(mLBP)
    tLBP = MathProgBase.getsolvetime(mLBP)
    sip_sto.LBP_time += tLBP
    sip_sto.LBD = LBDg
    sip_sto.xbar = xbar
    #  println("lower recall")
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
      return LBDg,UBDg,sip_sto
    else
      push!(P_LBD,pbar)
    end
    ##### upper bounding problem #####
    gUBP = x -> BndProb_reform(x,nothing,gSIP,P_UBD,eps_g)
    gL_UBP = [-Inf for i=1:length(P_UBD)]
    gU_UBP = [0.0 for i=1:length(P_UBD)]
    mUBP = deepcopy(MathProgBase.NonlinearModel(SIPopt.UBP_Opt))
    MathProgBase.loadproblem!(mUBP, nx, length(P_UBD), X_low, X_high,
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
