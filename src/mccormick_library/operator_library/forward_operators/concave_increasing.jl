for opMC in (:log, :log2, :log10, :log1p, :acosh, :sqrt)

   # get derivative of midcc for nonsmooth McCormick
   dop = diffrule(:Base, opMC, :midcc) # Replace with cv ruleset

   # creates expression for nonsmooth McCormick operator
   MC_exp = quote
              xIntv = ($opMC)(Intv(x))
              xL = lo(x)
              xU = hi(x)
              xLc = lo(xIntv)
              xUc = hi(xIntv)
              midcc,cc_id = mid3(cc(x),cv(x),xU)
              midcv,cv_id = mid3(cc(x),cv(x),xL)
              dcv = (xUc > xLc) ? (xUc-xLc)/(xU-xL) : 0.0
              convex = dcv*(midcv-xL)+ xLc
              concave = ($opMC)(midcc)
              concave_grad = mid_grad(cc_grad(x), cv_grad(x), cc_id)*$dop
              convex_grad = mid_grad(cc_grad(x), cv_grad(x), cv_id)*dcv
              convex, concave, convex_grad, concave_grad = cut(xLc,xUc,convex,concave,convex_grad,concave_grad)
            end

    # creates expression for smooth McCormick operator
    dop = diffrule(:Base, opMC, :(cv(x)))
    dMC_exp = quote
               xIntv = ($opMC)(Intv(x))
               xL = lo(x)
               xU = hi(x)
               xLc = lo(xIntv)
               xUc = hi(xIntv)
               midcc,cc_id = mid3(cc(x),cv(x),xU)
               midcv,cv_id = mid3(cc(x),cv(x),xL)
               dcv = (xUc > xLc) ? (hi(xIntv)-lo(xIntv))/(xU-xL) : 0.0
               if (x2-x1) == 0.0
                  convex = xLc
               else
                  convex = xLc*((xU-midcv)/(xU-xL)) + xUc*((midcv-xL)/(xU-xL))
               end
               concave = ($opMC)(midcc)
               convex_grad = ((hi(xIntv)-lo(xIntv))/(xU-xL))*cv_grad(x)
               concave_grad = ($dop)*cc_grad(x)
              end

    # calculates cse optimized expressions
    cse_MC_exp = MC_exp
    cse_diffMC_exp = dMC_exp

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
