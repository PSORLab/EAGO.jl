for opMC in (:exp, :exp2, :exp10, :expm1)
   opMC_kernel = Symbol(String(opMC)*"_kernel")
   dop = diffrule(:Base, opMC, :midcv) # Replace with cv ruleset
   MCexp = quote
              xL = x.Intv.lo
              xU = x.Intv.hi
              xLc = z.lo
              xUc = z.hi
              midcc, cc_id = mid3(x.cc, x.cv, xU)
              midcv, cv_id = mid3(x.cc, x.cv, xL)
              concave = xUc
              (xUc > xLc) && (concave = xLc*((xU - midcc)/(xU - xL)) + xUc*((midcc - xL)/(xU - xL)))
              convex = ($opMC)(midcv)
              concave_grad = mid_grad(x.cc_grad, x.cv_grad, cc_id)*(xUc - xLc)/(xU - xL)
              convex_grad = mid_grad(x.cc_grad, x.cv_grad, cv_id)*$dop
              convex, concave, convex_grad, concave_grad = cut(xLc, xUc, convex, concave, convex_grad, concave_grad)
              return MC{N, NS}(convex, concave, z, convex_grad, concave_grad, x.cnst)
            end
    dop = diffrule(:Base, opMC, :(x.cv))
    dMCexp = quote
               xL = x.Intv.lo
               xU = x.Intv.hi
               xLc = z.lo
               xUc = z.hi
               midcc,cc_id = mid3(x.cc, x.cv, xU)
               midcv,cv_id = mid3(x.cc, x.cv, xL)
               concave = xUc
               (xUc > xLc) && (concave = xLc*((xU-midcc)/(xU-xL)) + xUc*((midcc-xL)/(xU-xL)))
               convex = ($opMC)(midcc)
               convex_grad = ($dop)*x.cv_grad
               concave_grad = ((xUc - xLc)/(xU - xL))*x.cc_grad
               return MC{N, Diff}(convex, concave, z, convex_grad, concave_grad, x.cnst)
              end

      @eval @inline ($opMC_kernel)(x::MC{N, NS}, z::Interval{Float64}) where {N} = $MCexp
      @eval @inline ($opMC_kernel)(x::MC{N, Diff}, z::Interval{Float64}) where {N} = $dMCexp
      @eval @inline ($opMC)(x::MC) = ($opMC_kernel)(x, ($opMC)(x.Intv))
end
