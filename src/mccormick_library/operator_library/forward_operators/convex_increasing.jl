for opMC in (:exp,:exp2,:exp10,:expm1)

   # get derivative of midcv for nonsmooth McCormick
   dop = DiffRules.diffrule(:Base, opMC, :midcv) # Replace with cv ruleset

   monotone = :increasing
   if monotone == :increasing
         eps_min = :(lo(xIntv))
         eps_max = :(hi(xIntv))
    end
   # creates expression for nonsmooth McCormick operator
   MC_exp = quote
              xIntv = ($opMC)(Intv(x))
              xL = lo(x)
              xU = hi(x)
              xLc = lo(xIntv)
              xUc = hi(xIntv)
              midcc,cc_id = mid3(cc(x),cv(x),xU)
              midcv,cv_id = mid3(cc(x),cv(x),xL)
              concave = xUc
              (xUc > xLc) && (concave = xLc*((xU-midcc)/(xU-xL)) + xUc*((midcc-xL)/(xU-xL)))
              convex = ($opMC)(midcv)
              concave_grad = mid_grad(cc_grad(x), cv_grad(x), cc_id)*(hi(xIntv)-lo(xIntv))/(xU-xL)
              convex_grad = mid_grad(cc_grad(x), cv_grad(x), cv_id)*$dop
              convex, concave, convex_grad, concave_grad = cut(xLc,xUc,convex,concave,convex_grad,concave_grad)
            end

    # creates expression for nonsmooth McCormick operator
    dop = DiffRules.diffrule(:Base, opMC, :(cv(x)))
    dMC_exp = quote
               xIntv = ($opMC)(Intv(x))
               xL = lo(x)
               xU = hi(x)
               xLc = lo(xIntv)
               xUc = hi(xIntv)
               midcc,cc_id = mid3(cc(x),cv(x),hi(xIntv))
               midcv,cv_id = mid3(cc(x),cv(x),lo(xIntv))
               concave = xUc
               (xUc > xLc) && (concave = xLc*((xU-midcc)/(xU-xL)) + xUc*((midcc-xL)/(xU-xL)))
               convex = ($opMC)(midcc)
               convex_grad = ($dop)*cv_grad(x)
               concave_grad = ((hi(xIntv)-lo(xIntv))/(xU-xL))*cc_grad(x)
              end

    # calculates cse optimized expressions
    cse_MC_exp = cse(MC_exp)
    cse_diffMC_exp = cse(dMC_exp)

    # combines expression into single McCormick operator & overloads method
    comb_MC = quote
                 if (MC_param.mu >= 1)
                    $cse_diffMC_exp
                 else
                    $cse_MC_exp
                 end
                 return MC{N}(convex, concave, xIntv, convex_grad, concave_grad, x.cnst)
              end

     # generates the function
     @eval ($opMC)(x::MC{N}) where N = $comb_MC
end
