EpsMin_DecrToIncr = :((xU < sw) ? xU : ((sw < xL) ? xL : sw))
EpsMin_IncrToDecr = :((fxU < fxL) ? xU : xL)
EpsMin = Dict{Symbol,Union{Symbol,Expr}}(:Increasing => :xL,
                           :Decreasing => :xU,
                           :DecrToIncr => EpsMin_DecrToIncr,
                           :IncrToDecr => EpsMin_IncrToDecr)

EpsMax_DecrToIncr = :((fxU > fxL) ? xU : xL)
EpsMax_IncrToDecr = :((xU < sw) ? xU : ((xL > sw) ? xL : sw))
EpsMax = Dict{Symbol,Union{Symbol,Expr}}(:Increasing => :xU,
                           :Decreasing => :xL,
                           :DecrToIncr => EpsMax_DecrToIncr,
                           :IncrToDecr => EpsMax_IncrToDecr)

```
    DefineRootExpr(md::Symbol,op::Symbol)

Generates the rhs of functions to used to find the anchor point of the McCormick relaxation.
# convex relaxation for convexoconcave (or concave relaxation for concavoconvex)  f(p) = f'(p)*(xU-p) - f(xU) + f(p) = 0 (eqn)
                                                                                  f'(p) = f"(p)*(xU - p) = 0 (d eqn)
# convex relaxation for convexoconcave (or concave relaxation for concavoconvex)  f(p) = f'(p)*(p-xL) - f(p) + f(xL) = 0  (eqn)
                                                                                  f'(p) = f"(p)*(p - xL) = 0 (d eqn)
```
function DefineRootExpr(md::Symbol,op::Symbol)
    dop = DiffRules.diffrule(md, op, :x)
    root = :($dop*(y-x) - ($op)(y) + ($op)(x))
    droot = Calculus.differentiate(root,:x)
    String_Env = String(md)*"_"*String(op)*"_env"
    String_DEnv = String(md)*"_"*String(op)*"_denv"
    Sym_Env = Symbol(String_Env)
    Sym_DEnv = Symbol(String_DEnv)
    root, droot, Sym_Env, Sym_DEnv
    rootdef = quote
                  @inline ($Sym_Env)(x::Float64,y::Float64,z::Float64) = $root
                  @inline ($Sym_DEnv)(x::Float64,y::Float64,z::Float64) = $droot
               end
    return rootdef
end

function GetOuterRelax(mod, cvx, op, midcc, midcv)
         diffcv = DiffRules.diffrule(mod, op, midcv)
         diffcc = DiffRules.diffrule(mod, op, midcc)
         if (cvx == :Convex || cvx == :Affine) # Pretty Much Done
            gen_dconcave = :((xU > xL) ? (fxU - fxL)/(xU - xL) : 0.0)
            gen_concave = :((xU > xL) ? ((fxU - fxL)/(xU - xL))*($midcc-xL) + fxL, 0.0  : xU, 0.0)
            gen_dconvex = diffcv
            gen_convex = :(($op)($midcv), 0.0)
         elseif (cvx == :Concave) # Pretty Much Done
            gen_dconcave = diffcc
            gen_concave = :(($op)(midcc), 0)
            gen_dconvex = :((xU > xL) ? (fxU - fxL)/(xU - xL) : 0.0)
            gen_convex = :((xU > xL) ? ((fxU - fxL)/(xU - xL))*($midcc-xL) + fxL , 0.0  : xL, 0.0)
         elseif (cvx == :Convexoconcave)
            dccMIX = :(($midcc < p_cc) ? ((xU > xL) ? (fxU - fxL)/(xU - xL) : 0.0) : $diffcc)
            dcvMIX = :(($midcv < p_cv) ? $diffcv : ((xU > xL) ? (fxU - fxL)/(xU - xL) : 0.0))
            ccMIX =
            cvMIX =
            gen_dconcave = :(xU < sw ? :((xU > xL) ? (fxU - fxL)/(xU - xL) : 0.0) : ((xL > sw) ? $diffcc: $dccMIX))
            gen_concave = :(xU < sw ? :((xU > xL) ? ((fxU - fxL)/(xU - xL))*($midcc - xL) + fxL , 0.0 : xU , 0.0) : ((xL > sw) ? :($op)(midcc), 0.0 : ($ccMIX)(midcc,xL,xU))
            gen_dconvex = :(xU < sw ? $diffcv : ((xL > sw) ? :((xU > xL) ? (fxU - fxL)/(xU - xL) : 0.0) : $dcvMIX))
            gen_convex = :(xU < sw ? :($op)($midcv), 0.0 : ((xL > sw) ? :((xU > xL) ? ((fxU - fxL)/(xU - xL))*($midcv-xL) + fxL, 0.0  : xU, 0.0) : ($cvMIX)(midcc,xL,xU)))
         elseif (cvx == :Concavoconvex)
            dccMIX = :(($midcc < p_cc) ? $diffcc : ((xU > xL) ? (fxU - fxL)/(xU - xL) : 0.0))
            dcvMIX = :(($midcv < p_cv) ? ((xU > xL) ? (fxU - fxL)/(xU - xL) : 0.0) : $diffcv)
            ccMIX =
            cvMIX =
            gen_dconcave = :(xU < sw ?  $diffcc : ((xL > sw) ? :((xU > xL) ? (fxU - fxL)/(xU - xL) : 0.0) : $dccMIX))
            gen_concave = :(xU < sw ? :($op)(midcc) : ((xL > sw) ? :((xU > xL) ? ((fxU - fxL)/(xU - xL))*($midcc - xL) + fxL : xU) : ($ccMIX)(midcc,xL,xU)))
            gen_dconvex = :(xU < sw ? :((xU > xL) ? (fxU - fxL)/(xU - xL) : 0.0) : ((xL > sw) ? $diffcv : $dcvMIX))
            gen_convex = :(xU < sw ? :($op)($midcv) : ((xL > sw) ? :((xU > xL) ? ((fxU - fxL)/(xU - xL))*($midcv-xL) + fxL  : xU) : $(cvMIX)(midcc,xL,xU)))
         else
            error("Convexity rules not defined. Cannot generate McCormick relaxations")
         end
         return gen_dconcave, gen_concave, gen_dconvex, gen_convex
end

function McCormickExpr(md, op, ar)
   quote
   # Generates rootfinding functions


   # get derivative of midcv for nonsmooth McCormick
   cvx, mono, sw, low, upp = McCormick.ConvexityRules.CV_TRAITS[$md,$op,$ar]
   dcv = DiffRules.diffrule($md, $op, :midcv)
   dcc = DiffRules.diffrule($md, $op, :midcc)
   println("ran to me 1")
   gen_eps_min = EpsMin[$op]
   gen_eps_max = EpsMax[$op]
   println("ran to me 2")
   MCterms = GetOuterRelax($mod, cvx, $op, midcc, midcv)

   #MC_exp =
   #dMC_exp =

   # calculates cse optimized expressions for nonsmooth and smooth McCormick operators
   #cse_MC_exp = cse(MC_exp)
   #cse_diffMC_exp = cse(dMC_exp)

   # combines expression into single McCormick operator & overloads method
   comb_MC = quote
                 xIntv = Intv(x)
                 xL, xU = lo(xIntv), hi(xIntv)
                 fxIntv = ($op)(xIntv)
                 fLx, fUx = lo(fxIntv), hi(fxIntv)
                 fxL, fxU = ($op)(xL), $(op)(xU)
                 eps_min = $gen_min
                 eps_max = $gen_max
                 midcc, cc_id = mid3(cc(x), cv(x), eps_max)
                 midcv, cv_id = mid3(cc(x), cv(x), eps_min)
                 if (MC_param.mu >= 1)
                   # $cse_diffMC_exp
                 else
                  #  $cse_MC_exp
                 end
                 return MC{N}(convex, concave, xIntv, convex_grad, concave_grad, x.cnst)
              end

     # generates the function
     @eval ($op)(x::MC{N}) where N = $comb_MC
  end
end
#=
quote
     concave =
     convex = ($op)(midcv)
     concave_grad = mid_grad(cc_grad(x), cv_grad(x), cc_id)*($gen_dconcave)
     convex_grad = mid_grad(cc_grad(x), cv_grad(x), cv_id)*($gen_dconvex)
     convex, concave, convex_grad, concave_grad = cut(xLc,xUc,convex,concave,convex_grad,concave_grad)
end
=#

@GenerateMcCormick(:exp)
#=
q1 = quote dop = DiffRules.diffrule(:Base, $opMC, :midcv) end
q2 = quote mono = cvtrait() end
q3 = quote convex = 1 end
q4 = quote smth = 1 end
println("$q1")
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
           midcc,cc_id = mid3(cc(x),cv(x),$eps_max)
           midcv,cv_id = mid3(cc(x),cv(x),$eps_min)
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
=#
