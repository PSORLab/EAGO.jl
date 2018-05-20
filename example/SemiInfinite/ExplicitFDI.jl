using IntervalArithmetic
using EAGO

function FDI_invApprox(x::Interval)
  Ap = 0.00005
  g = 9.8
  low = -(x^2)/(2*g*Ap^2)
  high = (x^2)/(2*g*Ap^2)
  if x.lo>0
    return high
  elseif x.hi<0
    return low
  else
    return Interval(0,max(low.hi,high.hi))
  end
end

function FDI_invApprox(x::Float64)
  Ap = 0.00005
  g = 9.8
  if x>=0
    return (x^2)/(2*g*Ap^2)
  else x<0
    return -(x^2)/(2*g*Ap^2)
  end
end

function xbnds(u,p)

         # parameters
         Ap = 0.00005
         g = 9.8
         z = []
         # bound z via McCormick operators
         z1 = (((u[1]/p[3])/Ap)^2)/(2.0*g)
         z3 = (1.0/(2.0*g))*((u[1]+u[2])/(p[4]*Ap))^2
         z2 = (-((u[2]-p[4]*Ap*sqrt(2*g*z3))/(p[5])/Ap)^2)/(2.0*g)
        # z4 = (((u[1]*p[1]/p[3])/Ap)^2)/(2.0*g)
        # z6 = (1.0/(2.0*g))*((p[1]*u[1]+u[2])/(p[4]*Ap))^2
        # z5 = (-((u[2]-p[4]*Ap*sqrt(2*g*z6))/(p[5])/Ap)^2)/(2.0*g)
      #   z7 = (((u[1]/p[3])/Ap)^2)/(2.0*g)
      #   z9 = (1.0/(2.0*g))*((u[1]+u[2])/(p[4]*(Ap+3.14159*p[2]^2)))^2
      #   z8 = u[2]-p[4]*(Ap+3.14159*p[2]^2)*sqrt(2.0*g*z9)/(p[5])
        # push!(z,z1,z2,z3,z4,z5,z6,z7,z8,z9)
         push!(z,z1,z2,z3)
         # rotate coordinates back to X & intersects

         x = []
         x1 = min(max(z1-z2+z3,0.0),0.75)
         x2 = min(max(z3,0.0),0.75)
         x3 = min(max(z3-z2,0.0),0.75)
        # x4 = min(max(z4-z5+z6,0.0),0.75)
         #x5 = min(max(z6,0.0),0.75)
        # x6 = min(max(z6-z5,0.0),0.75)
         #x7 = min(max(z7-z8+z9,0.0),0.75)
         #x8 = min(max(z9,0.0),0.75)
        # x9 = min(max(z9-z8,0.0),0.75)
         #x1 = z1-z2+z3
         #x2 = z3
         #x3 = z3-z2
         #x4 = z4-z5+z6
         #x5 = z6
         #x6 = z6-z5
         #x7 = z7-z8+z9
        # x8 = z9
        # x9 = z9-z8
         #push!(x,x1,x2,x3,x4,x5,x6,x7,x8,x9)
         push!(x,x1,x2,x3)#,x4,x5,x6,x7,x8,x9)
         return x
end

function SIP_g(un,p)
  u = un[1:2]
  xbnd = xbnds(u,p)
  r_xbnds = [zero(xbnd[1]) for i=1:3]
  #r_xbnds[:,1] = xbnd[1:3]
  #r_xbnds[:,2] = xbnd[4:6]
  #r_xbnds[:,3] = xbnd[7:9]
  out = un[3]
  for i=1:3
        out = out - (r_xbnds[i])^2
  end
  return out
end

# sets up SIP objective
function SIP_f(un)
  return -un[3]
end

#ex_SIP_X = [(1E-5)..(1E-4),(1E-5)..(1E-4),(-10.0)..(10.0)]
ex_SIP_X = [(1E-4)..(1.00001E-4),(1E-4)..(1.00001E-4),(-10.0)..(10.0)]
ex_SIP_P = [(0.54)..(0.66),(0.0005)..(0.005),(0.85)..(1.15),(0.65)..(0.95),(0.85)..(1.15)]

SIPopt = SIP_opts()
sep1lu = EAGO_NLPSolver(LBD_func_relax = "NS-STD-Off",
                        LBDsolvertype = "LP",
                        probe_depth = -1,
                        variable_depth = 1000,
                        DAG_depth = -1,
                        STD_RR_depth = -1,
                        UBDsolvertype= "Ipopt",
                        validated = true,
                        atol = 1E-4,
                        rtol = 1E-4)
sep1lu.BnBSolver.Verbosity = "None"
sep1in = EAGO_NLPSolver(LBD_func_relax = "NS-STD-Off",
                        LBDsolvertype = "LP",
                        probe_depth = -1,
                        variable_depth = 1000,
                        DAG_depth = -1,
                        STD_RR_depth = -1,
                        UBDsolvertype= "Ipopt",
                        validated = true,
                        atol = 1E-4,
                        rtol = 1E-4)
sep1in.BnBSolver.Verbosity = "None"
SIPopt.LLP_Opt = sep1in
SIPopt.LBP_Opt = sep1lu
SIPopt.UBP_Opt = sep1lu
ex_out = Explicit_SIP_Solve(SIP_f,SIP_g,ex_SIP_X,ex_SIP_P,SIPopt)
