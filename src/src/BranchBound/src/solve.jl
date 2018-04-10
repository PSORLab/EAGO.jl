"""
    solveBnB!(x::BnBSolver,y::BnBModel)

Solves the branch and bound problem with the input model and solver object.
"""
function solveBnB!(x::BnBSolver,y::BnBModel)

  # initializes counts
  k_int::Int64 = 0
  k_nod::Int64 = length(y.LBD)

  # terminates when max nodes or iteration is reach, or when node stack is empty
  while (x.Term_Check(x,y,k_int))

    # fathom nodes with lower bound greater than global upper bound
    y.LBDg = minimum(y.LBD)
    push!(y.LBDg_hist,minimum(y.LBD))

    # selects node (deletion included)
    nsBox,LBDn,UBDn,id,pos = x.Node_Select(y)

    # prints node in full verbosity mode
    print_node!(x,id,LBDn,nsBox)

    # solves preprocessing/LBD/UBD/postprocessing once to get timing right
    feas_Pre::Bool = true
    feas_Post::Bool = true
    if (k_int == 0)
      println("pre-check")
      nsBox1 = copy(nsBox)
      yUBDg1 = copy(y.UBDg)
      yLBDg1 = copy(y.LBDg)
      kint1 = copy(k_int)
      pos1 = copy(pos)
      xopt1 = copy(x.opt)
      x.Preprocess(feas_Pre,nsBox1,yUBDg1,kint1,pos1,xopt1,LBDn,UBDn,x,y)
      println("preprocess-check")
      D_valt, LBD_solt, LBD_feast, temp_objtL = x.Lower_Prob(nsBox1,kint1,pos1,xopt1,yUBDg1)
      println("lower-check")
      UBD_valt, UBD_solt, UBD_feast, temp_objtU = x.Upper_Prob(nsBox1,kint1,pos1,xopt1,yUBDg1)
      println("upper-check")
      x.Postprocess(feas_Post,nsBox1,kint1,pos1,xopt1,
                    temp_objtL,temp_objtU,yLBDg1,yUBDg1)
      println("postprocess-check")
    end
    feas_Pre = true
    feas_Post = true

    # performs prepocessing and times
    tic()
    feas_Pre,nsBox = x.Preprocess(feas_Pre,nsBox,y.UBDg,k_int,pos,x.opt,LBDn,UBDn,x,y)
    push!(y.Pretime,y.Pretime[end]+toq())

    UBD_feas = false
    if (feas_Pre)
      # solves & times lower bounding problem
      tic()
      LBD_val,LBD_sol,LBD_feas,temp_objL = x.Lower_Prob(nsBox,k_int,pos,x.opt,y.UBDg)
      push!(y.LBDgtime,y.LBDgtime[end]+toq())
      y.lbcnt += 1
      print_results!(x,LBD_val,LBD_sol,LBD_feas,true)
      int_info = LBD_val,LBD_sol,LBD_feas,temp_objL,nsBox

      # checks for infeasibility stores solution
      if (LBD_feas)
      #  println("ran LBD feas")
      #  println("boolean: $(x.converged(x,y.UBDg,LBD_val))")
        if (~x.converged(x,y.UBDg,LBD_val))

          # solves & times upper bounding problem
        #  println("ran Upper bound")
          tic()
          UBD_val,UBD_sol,UBD_feas,temp_objU = x.Upper_Prob(nsBox,k_int,pos,x.opt,y.UBDg)
          push!(y.UBDgtime,y.UBDgtime[end]+toq())
          y.ubcnt += 1
          print_results!(x,UBD_val,UBD_sol,UBD_feas,false)

          # fathoms by value dominance
          if (UBD_feas)
            # update to handle equality constraints
            if (y.lastgap>(UBD_val-LBD_val) && (UBD_val - y.UBDg < x.BnB_atol/10.0))
              y.lastgap = UBD_val-LBD_val
              y.soln = UBD_sol
            end
            if (UBD_val < y.UBDg)
              y.feas_fnd = true
              y.first_num = y.lbcnt
              y.UBDg = UBD_val
              y.soln = UBD_sol
              push!(y.UBDg_hist,UBD_val)
            else
              push!(y.UBDg_hist,y.UBDg)
            end
          else
            push!(y.UBDg_hist,y.UBDg_hist[end])
          end

          # performs post processing and times
          tic()
          feas_Post,nsBox_Post = x.Postprocess(feas_Post,nsBox,k_int,pos,x.opt,temp_objL,temp_objU,y.LBDg,y.UBDg)
          push!(y.Posttime,y.Posttime[end]+toq())

          # branch if criteria met
          if (feas_Post)
            if x.Repeat_Chk(x,y,nsBox,nsBox_Post)
              BM_Single!(x,y,LBD_val,UBD_val,nsBox,pos)
            else
              # branch the nodes & stores nodes to stack
              Y1,Y2 = x.Bisect_Func(x,y,nsBox)
              x.Branch_Sto(x,y,LBD_val,UBD_val,Y1,Y2,pos)
            end
          end

        elseif (~x.exhaust && y.feas_fnd)
        end# end converged check
      end # end LBD feasibility check

      # prints relative statistics for the iteration
      fathom!(y)
    else
      LBD_val = -Inf
      LBD_feas = false
      UBD_feas = false
    end

    print_int!(x,k_int,length(y.LBD),id,LBD_val,y.LBDg,y.UBDg,LBD_feas,UBD_feas)

    k_int+=1

  end
  y.soln_val = y.UBDg

  # prints the solution
  if ((x.Verbosity == "Full")||(x.Verbosity == "Normal"))
    print_sol!(x,y,y.ubcnt,y.lbcnt,y.UBDgtime[end],y.LBDgtime[end])
  end
end
