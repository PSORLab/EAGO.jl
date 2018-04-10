# Current version has flaws since update to methods will fixed in future
#=
"""
--------------------------------------------------------------------------------
Function: Generalized_Param_Bisection
--------------------------------------------------------------------------------
Description:
Uses the interval arithmetic presented in ValidatedNumerics to bound unique
implicit functions x:P->X defined by a system of equations h(z,p)=0 via
Gauss-Siedel Newton.
--------------------------------------------------------------------------------
Inputs:
* Xin:       Vector{Interval{Float64}} - Bounds for dependent variables
* Pin:       Vector{Interval{Float64}} - Bounds for independent variables
* h:         function - Equations defining potential implicit function
* hj:        function - Jacobian of h(x,p) with respect to x
* opt:       Param_Bisect_Opts - Parametric Bisection Options
--------------------------------------------------------------------------------
Returns:
A tuple (sol_set, ind_set, NIT, hist) where
* sol_set:     Array - Interval boxes known to contain unique implicit function
* ind_set:     Array - Interval boxes known to contain unique implicit function
* NIT:         Int64 - Final iteration number
* hist:        Array - Storage object with info on problem solution history
--------------------------------------------------------------------------------
"""
#=
function Generalized_Param_Bisection(Xin,Pin,h,hj,opt)

  # unpacks options for general bisection
  DAGflag::Bool = opt.DAGflag
  LPflag::Bool = opt.LPflag
  max_iter = opt.kmax_main
  max_iter_cntr = opt.kmax_cntr
  style = opt.style
  display = opt.display
  ptol = opt.ptol
  etol = opt.etol
  rtol = opt.rtol
  pbisect = opt.p_rel_bisect
  Pstart = copy(Pin)
  bp = 1

  if (DAGflag)
    println("started generating DAG parameters")
    DAGr = opt.DAGpass
    DAGh = opt.DAGh
    DAGg = opt.DAGg
    DAGsym = opt.DAGsym
    DAGpack = vcat(Xin,Pin)
    DAGparam = Contractor_Params(DAGh,DAGg,DAGpack,DAGsym)
    println("generated DAG parameters")
  end

  inclusionHigh = [false for i=1:length(Xin)]
  inclusionLow = [false for i=1:length(Xin)]
  opt = [max_iter_cntr,etol,rtol,style]

  # sets dimensions
  np = length(Pin)
  nx = length(Xin)

  # initializes count variables
  NIT = 0
  NMEX = 0
  NPINMEX = 0
  NS = 0

  # creates working stack & storage objects, pushes initial box onto stack
  sol_set = []
  ind_set = []
  stack = []
  hist = []
  push!(stack,[Xin,Pin])
  disp = 0
  # calculates DAG contractor parameters

  # Check Termination
  while (~isempty(stack) && (NIT < max_iter))
    #if (display)
      println("Iteration  ||  Remaining  ||  Solution Size  ||  Indeterminat Size")
      println(NIT, "  ||  ",length(stack),"  ||  ",length(sol_set),"  ||  ",length(ind_set))
    #end
    #### Resets the tests variables and unpacks the stack ####
    Iflag = false
    Eflag = false
    PIflag = false
    PIcert = false
    eDflag = false
    cnode = pop!(stack)
    Xw = copy(cnode[1])
    X0 = copy(cnode[1])
    Pw = copy(cnode[2])
    println("Current Node (X,P): ", Xw, " ",Pw)

    #### Check for exclusion via Miranda's Test ####
    Eflag = MirandaExc(h,Xw,Pw,Eflag)
    if (Eflag)
      println("Miranda Exclusion:  ",Eflag,"  h(Z,P):  ",h(Xw,Pw))
    end
    Eflag && (NMEX += 1)
    disp = 1
    histout = []

    if (~Eflag)
      #=
      if (DAGflag)
        DAGpack = vcat(Xw,Pw)
        DAGpack = Contractor_Eval(DAGh,DAGg,DAGpack,DAGsym,DAGr,DAGparam)
        #println("typeof(DAGpack): ",typeof(DAGpack))
        #println("typeof(Xw): ",typeof(Xw))
        for i=1:length(DAGpack)
          if (isempty(DAGpack[i]))
            println("Fathom by DAG")
            Eflag = true
            break
          end
        end
        for i=1:length(Xw)
          Xw[i] = DAGpack[i]
        end
        count = 1
        for i=(length(Xw)+1):(length(Xw)+length(Pw))
          Pw[count] = DAGpack[i]
          count += 1
        end
      end

      println("ran DAG contractor")
    end

    if (~Eflag)
      #println("parametric iterations start")
      if (style == "NewtonGS")
        Xw,Xw2,Eflag,Iflag,eDflag,inclusionLow,inclusionHigh = PI_NewtonGS(Xw,Pw,hj,h,opt,Eflag,Iflag,eDflag)
        #println("Xw:   ", Xw)
        #println("Xw2:   ", Xw2)
        #println("Eflag:   ", Eflag)
        #println("Iflag:   ", Iflag)
        #println("eDflag:   ", eDflag)
        #println("inclusionLow:   ", inclusionLow)
        #println("inclusionHigh:   ", inclusionHigh)
        if eDflag
          push!(stack,[Xw2,Pw])
        end
      end
      if (style == "KrawczykCW")
      #  println("Xw Initial: ",Xw)
        Xw,Eflag,Iflag,inclusionLow,inclusionHigh =  PI_KrawczykCW(Xw,Pw,hj,h,opt,Eflag,Iflag)
      #  println("Xw Final: ",Xw)
      end
      #println("eDflag: ", eDflag)
    #  println("Eflag: ", Eflag)
      #println("Iflag: ", Iflag)
    #  println("PIflag: ", PIflag)
    #  println("PIcert: ", PIcert)

      #println("~Iflag: ", ~Iflag)
    #  println("~PIflag: ", ~PIflag)
      #println("(~PIflag || PIcert): ", (~PIflag || PIcert))
    #  println("(~Iflag && (~PIflag || PIcert)): ", (~Iflag && (~PIflag || PIcert)))

      println("start boundary test:")
      if (Eflag==true)
        disp = 2
        NPINMEX += 1
      elseif (~Iflag && (~PIflag || PIcert))
        PIcert,PIflag,Iflag,Eflag = BoundaryTest(h,hj,X0,Xw,Pw,opt,PIcert,PIflag,
                                                 Iflag,Eflag,inclusionLow,inclusionHigh)
        println("performed boundary test")
        println("boundary test results")
        println("Eflag: ", Eflag)
        println("Iflag: ", Iflag)
        println("PIflag: ", PIflag)
        println("PIcert: ", PIcert)
      end

      if (Iflag)
        println("start re-refine")
        if (style == "NewtonGS")
          #println("trace 4")
          Xw,Xw2,Eflag,Iflag,eDflag,inclusionLow,inclusionHigh =  MC_NewtonGS(Xw,Pw,hj,h,opt,Eflag,Iflag,eDflag)
          #println("trace 4a")
          #println("Xw:   ", Xw)
          #println("Xw2:   ", Xw2)
          #println("Eflag:   ", Eflag)
          #println("Iflag:   ", Iflag)
          #println("eDflag:   ", eDflag)
          #println("inclusionLow:   ", inclusionLow)
          #println("inclusionHigh:   ", inclusionHigh)
        end
        if (style == "KrawczykCW")
          #println("trace 4")
          Xw,Eflag,Iflag,inclusionLow,inclusionHigh =  MC_KrawczykCW(Xw,Pw,hj,h,opt,Eflag,Iflag)
          #println("trace 4a")
        end
        disp = 3
        push!(sol_set,[Xw,Pw])
        NS += 1
      elseif (~Eflag)
        println("start XP bisection")
        NextNode = copy(Xw)
        NextNodeP =  copy(Pw)
        Ptemp =  copy(Pw)
        PIflag = false
        bp,CNodeX,CNodeP,NNodeX,NNodeP,Ptemp = XP_Bisection(h,hj,Xw,Pw,NextNode,NextNodeP,Ptemp,
                                                            PIflag,Iflag,Eflag,ptol,pbisect,Pstart)
        println("trace 5 - Bisection Results")
        println("bp:  ", bp)
        println("CNodeX:  ", CNodeX)
        println("CNodeP:  ", CNodeP)
        println("NNodeX:  ", NNodeX)
        println("NNodeP:  ", NNodeP)
        println("Ptemp:  ", Ptemp)
        if (bp == 2)
          disp = 4
          push!(stack,[CNodeX,CNodeP],[CNodeX,NNodeP])
        elseif (bp == 3)
          disp = 5
          push!(stack,[CNodeX,CNodeP],[CNodeX,NNodeP],[CNodeX,Ptemp])
        elseif (bp == 1)
          disp = 6
          push!(stack,[CNodeX,CNodeP],[NNodeX,CNodeP])
        else
          disp = 7
          push!(ind_set,[CNodeX,CNodeP])
        end
      end
    end
    #println("trace end loop")
    push!(hist,[NIT,cnode,disp,histout])
    NIT += 1
  end
  #println("trace end")
  return sol_set, ind_set, NIT,hist
end
=#
