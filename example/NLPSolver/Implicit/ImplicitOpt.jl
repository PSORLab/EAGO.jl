workspace()
using EAGO
using JuMP
using MathProgBase
#using IntervalArithmetic
#using ForwardDiff
using Ipopt
#using NLopt

# set Lowerbound type
LBD_fr = "NS-STD"
LBD_ps = "LP"
# set Interval Contractor type
intv_style = "KrawczykCW"
#intv_flag =
# set MC Contractor type
MC_style = "KrawczykCW"
#MC_flag =

###############################################################################
# Test Problem #1 (Solves Implicit Example 3 with LP Lower Bounds, In Place)
# Expected Minima: -2.61E+06
###############################################################################

opt1 = Any[5    #Int64: Number of iterations
       1.0E-6 #Float64: Tolerance for equality of
       1.0E-6 #Float64: Add Interval(1E-8,1E8) to add to M[i,i] when
              #         processing extended interval division.
      ]
#=
# loads the implicit solver options
f1(x,p) = x[1]^2 - (p[1]^2)*x[2] + p[2]
function h1!(hout::Vector{T},x,p) where T
       hout[1] = one(p[1])*((3.25-x[1])/p[1]-x[3])
       hout[2] = one(p[1])*((x[1]/p[2]-x[3]))
       hout[3] = one(p[1])*(x[2]-(x[1]^2)/(one(p[1])+x[1]^2))
end
function hj1!(hout::T,x,p) where T
       hout[1,1] = -one(p[1])/p[1]
       hout[1,2] = zero(p[1])
       hout[1,3] = -one(p[1])
       hout[2,1] = one(p[1])/p[2]
       hout[2,2] = zero(p[1])
       hout[2,3] = -one(p[1])
       hout[3,1] = one(p[1])*(-2.0*x[1]/(one(p[1])+x[1]^2)^2)
       hout[3,2] = one(p[1])
       hout[3,3] = zero(p[1])
end

function LBD_func(i)
          YL = [-30, -30, -30, 1800, 900]
          return YL[i]
end
function UBD_func(i)
          YU = [30, 30, 30, 2200, 1100]
          return YU[i]
end

function h2(x,p)
       return [(3.25-x[1])/p[1]-x[3]; x[1]/p[2]-x[3]; x[2]-(x[1]^2)/(one(p[1])+x[1]^2)]
end
function hj2(x,p)
       return [-one(p[1])/p[1] zero(p[1]) -one(p[1]);
       one(p[1])/p[2] zero(p[1]) -one(p[1]);
       -2.0*x[1]/(one(p[1])+x[1]^2)^2 one(p[1]) zero(p[1])]
end

imp_opt = ImplicitSolver()
imp_opt.opts.style = MC_style #"KrawczykCW" # sets MC contractor style
imp_opt.opts.kmax = 3
imp_opt.h = h2
imp_opt.hj = hj2
imp_opt.f = f1
imp_opt.g = [] # inherites bounds from jump model
imp_opt.nx = 3
imp_opt.flag = true
imp_opt.Intv_Cntr =  intv_style #"NewtonGS" # sets interval style
imp_opt.ParamInt = opt1
imp_opt.numConstr = 0
imp_opt.gL_Loc = []
imp_opt.gU_Loc = []
imp_opt.gL = []
imp_opt.gU = []

jumpmodel = Model(solver=EAGO_NLPSolver(ImplicitOpts = imp_opt,
                                        LBD_func_relax = LBD_fr,
                                        LBDsolvertype = LBD_ps,
                                        probe_depth = -1,
                                        variable_depth = -1,
                                        STD_RR_depth = -1,
                                        DAG_depth = -1))
x = @variable(jumpmodel, [i=1:5], lowerbound=LBD_func(i), upperbound=UBD_func(i))
@NLconstraint(jumpmodel, (3.25-x[1])/x[4]-x[3] == 0.0 )
@NLconstraint(jumpmodel, x[1]/x[5]-x[3] == 0.0 )
@NLconstraint(jumpmodel, x[2]-(x[1]^2)/(1.0+(x[1]^2)) == 0.0 )
@NLobjective(jumpmodel, Min, x[1]^2-(x[4]^2)*x[2]+x[5])
status = solve(jumpmodel)

#######################################################################################
# Test Problem #2 (Solves Implicit Example 3 with NLP Lower Bounds, w/Constr)
# Expected Minima: -2.07E+06
#######################################################################################

opt1 = Any[3    #Int64: Number of iterations
       1.0E-8 #Float64: Tolerance for equality of
       1.0E-8 #Float64: Add Interval(1E-8,1E8) to add to M[i,i] when
              #         processing extended interval division.
      ]

# loads the implicit solver options
f1(x,p) = x[1]^2 - (p[1]^2)*x[2] + p[2]
function h1!(hout,x,p)
       hout[1] = (3.25-x[1])/p[1]-x[3]
       hout[2] = x[1]/p[2]-x[3]
       hout[3] = x[2]-(x[1]^2)/(1+x[1]^2)
end
function hj1!(hout,x,p)
       hout[1,1] = -1.0/p[1]
       hout[1,2] = 0.0
       hout[1,3] = -1.0
       hout[2,1] = 1.0/p[2]
       hout[2,2] = 0.0
       hout[2,3] = -1.0
       hout[3,1] = -2.0*x[1]/(1.0+x[1]^2)^2
       hout[3,2] = 1.0
       hout[3,3] = 0.0
end

function LBD_func(i)
          YL = [-30, -30, -30, 1800, 900]
          return YL[i]
end
function UBD_func(i)
          YU = [30, 30, 30, 2200, 1100]
          return YU[i]
end

function h2(x,p)
       return [(3.25-x[1])/p[1]-x[3]; x[1]/p[2]-x[3]; x[2]-(x[1]^2)/(1+x[1]^2)]
end
function hj2(x,p)
       return [-1/p[1] 0 -1;
       1/p[2] 0 -1;
       -2*x[1]/(1+x[1]^2)^2 1 0]
end

function gmod2(x,p)
return [-(2*x[1]-x[3]-5),
          p[1]-1.7*p[2]]
end

imp_opt1 = ImplicitSolver()
imp_opt1.opts.style = MC_style #"KrawczykCW" # sets MC contractor style
imp_opt1.h = h2
imp_opt1.hj = hj2
imp_opt1.f = f1
imp_opt1.g = gmod2 # inherites bounds from jump model
imp_opt1.nx = 3
imp_opt1.flag = true
imp_opt1.Intv_Cntr = intv_style
imp_opt1.ParamInt = opt1
imp_opt1.numConstr = 0
imp_opt1.gL_Loc = []
imp_opt1.gU_Loc = [1 2]
imp_opt1.gL = [-Inf -Inf]
imp_opt1.gU = [0.0, 0.0]

jumpmodel2 =Model(solver=EAGO_NLPSolver( ImplicitOpts = imp_opt1,
                                         LBD_func_relax = LBD_fr,
                                         LBDsolvertype = LBD_ps,
                                         UBDsolvertype = "SNOPT",
                                         probe_depth = -1,
                                         variable_depth = -1,
                                         STD_RR_depth = -1,
                                         DAG_depth = -1))
x = @variable(jumpmodel2, [i=1:5], lowerbound=LBD_func(i), upperbound=UBD_func(i))
@NLconstraint(jumpmodel2, (3.25-x[1])/x[4]-x[3] == 0.0 )
@NLconstraint(jumpmodel2, x[1]/x[5]-x[3] == 0.0 )
@NLconstraint(jumpmodel2, x[2]-(x[1]^2)/(1+(x[1]^2)) == 0.0 )
@NLconstraint(jumpmodel2, 2*x[1]-x[3]-5 >= 0.0 )
@NLconstraint(jumpmodel2, x[4]-1.7*x[5] <= 0.0 )
@NLobjective(jumpmodel2, Min, x[1]^2-(x[4]^2)*x[2]+x[5])
status = solve(jumpmodel2)
=#

###############################################################################
# Test Problem #2 (Solves Quadratic, Constrained)
# Expected Minima: 20.25
###############################################################################
function h3(x,p)
       return [x[1]^2 + p[1]*x[1] + 4.0]
end

function hj3(x,p)
       return [2.0*x[1] + p[1]]
end

function f3(x,p)
       return (x[1]^2)*(p[1]^2)
end

function gmod3(x,p)
       return [9.0*x[1] + x[2]]
end

function LBD1_func(i)
          YL = [-0.78, 6.0]
          return YL[i]
end
function UBD1_func(i)
          YU = [-0.4, 9.0]
          return YU[i]
end

imp_opt1 = ImplicitSolver()
imp_opt1.h = h3
imp_opt1.hj = hj3
imp_opt1.f = f3
imp_opt1.g = gmod3 # inherites bounds from jump model
imp_opt1.nx = 1
imp_opt1.flag = true
imp_opt1.numConstr = 0
imp_opt1.gL_Loc = []
imp_opt1.gU_Loc = [1]
imp_opt1.gL = [-Inf]
imp_opt1.gU = [0.0]

jumpmodel3 = Model(solver=EAGO_NLPSolver( ImplicitOpts = imp_opt1,
                                         LBD_func_relax = LBD_fr,
                                         LBDsolvertype = LBD_ps,
                                         UBDsolvertype = "SNOPT",
                                         probe_depth = -1,
                                         variable_depth = -1,
                                         STD_RR_depth = -1,
                                         DAG_depth = -1))
x = @variable(jumpmodel3, [i=1:2], lowerbound=LBD1_func(i), upperbound=UBD1_func(i))
@NLconstraint(jumpmodel3, x[1]^2 + x[2]*x[1] + 4 == 0.0 )
@NLconstraint(jumpmodel3, 9*x[1] + x[2] <= 0.0 )
@NLobjective(jumpmodel3, Min, (x[1]^2)*(x[2]^2))
status = solve(jumpmodel3)


###############################################################################
# Test Problem #2 (Solves Quadratic, UnConstrained)
# Expected Minima: 17.8
###############################################################################
#=
function h3(x,p)
       return [x[1]^2 + p[1]*x[1] + 4.0]
end

function hj3(x,p)
       return [2*x[1] + p[1]]
end

function f3(x,p)
       return (x[1]^2)*(p[1]^2)
end

function gmod3(x,p)
       return [9*x[1] + x[2]]
end

function LBD1_func(i)
          YL = [-0.78, 6]
          return YL[i]
end
function UBD1_func(i)
          YU = [-0.4, 9]
          return YU[i]
end

imp_opt1 = ImplicitSolver()
imp_opt1.opts.style = MC_style #"KrawczykCW" # sets MC contractor style
imp_opt1.h = h3
imp_opt1.hj = hj3
imp_opt1.f = f3
imp_opt1.g = [] # inherites bounds from jump model
imp_opt1.nx = 1
imp_opt1.flag = true
imp_opt1.Intv_Cntr = intv_style
imp_opt1.ParamInt = opt1
imp_opt1.numConstr = 0
imp_opt1.gL_Loc = []
imp_opt1.gU_Loc = []
imp_opt1.gL = []
imp_opt1.gU = []

jumpmodel4 = Model(solver=EAGO_NLPSolver(ImplicitOpts = imp_opt1,
                                         LBD_func_relax = LBD_fr,
                                         LBD_problem_relax = LBD_pr,
                                         LBD_problem_solver = LBD_ps,
                                         UBD_func_relax = "Original",
                                         UBD_problem_relax = "NLP1",
                                         UBD_problem_solver = "SNOPT",
                                         probe_depth = -1,
                                         variable_depth = 1000,
                                         STD_RR_depth = -1,
                                         DAG_depth = -1))
x = @variable(jumpmodel4, [i=1:2], lowerbound=LBD1_func(i), upperbound=UBD1_func(i))
@NLconstraint(jumpmodel4, x[1]^2 + x[2]*x[1] + 4 == 0.0 )
@NLobjective(jumpmodel4, Min, (x[1]^2)*(x[2]^2))
status = solve(jumpmodel4)

###############################################################################
# Test Problem #5 (Reactor Separator)
# Expected Minima: 234295
###############################################################################
function hrs(z,p)
    # translate the inputs
    V = p[1]
    k1 = p[2]
    k2 = p[3]
    F1 = z[1]
    F2 = z[2]
    F3 = z[3]
    F4 = z[4]
    F6 = z[5]
    F7 = z[6]
    y3A = z[7]
    y3B = z[8]
    y3C = z[9]
    y4B = z[10]
    y4C = z[11]

    F5 = 50

    # parameters
    VA = 8.937E-2
    VB = 1.018E-1
    VC = 1.130E-1

    # intermediate calculations
    r1 = k1*y3A/(y3A*VA + y3B*VB + y3C*VC)
    r2 = k2*y3B/(y3A*VA + y3B*VB + y3C*VC)

    hout = zeros(typeof(p[1]),11)
    hout[1] = F1 + F7 - F2
    hout[2] = F2 - y3A*F3 - r1*V
    hout[3] = (r1+r2)*V - F5
    hout[4] = r2*V - y3C*F3
    hout[5] = 1.0 - y3A - y3B - y3C
    hout[6] = F3 - F4 - F7
    hout[7] = y3A*F3 - F7
    hout[8] = y3B*F3 - y4B*F4
    hout[9] = F4 - F5 - F6
    hout[10] = y4B*F4 - F5
    hout[11] = y4C*F4 - F6
    return hout
end

function hjrs(z,p)
    # translate the inputs
    V = p[1]
    k1 = p[2]
    k2 = p[3]
    F1 = z[1]
    F2 = z[2]
    F3 = z[3]
    F4 = z[4]
    F6 = z[5]
    F7 = z[6]
    y3A = z[7]
    y3B = z[8]
    y3C = z[9]
    y4B = z[10]
    y4C = z[11]

    # parameters
    VA = 8.937E-2
    VB = 1.018E-1
    VC = 1.130E-1

    # intermediate calculations
    r1 = k1*y3A/(y3A*VA + y3B*VB + y3C*VC)
    r2 = k2*y3B/(y3A*VA + y3B*VB + y3C*VC)

    r1y3A = k1*(y3B*VB + y3C*VC)/(y3A*VA + y3B*VB + y3C*VC)^2
    r1y3B = -k1*y3A*VB/(y3A*VA + y3B*VB + y3C*VC)^2
    r1y3C = -k1*y3A*VC/(y3A*VA + y3B*VB + y3C*VC)^2
    r2y3A = k2*(y3A*VA + y3C*VC)/(y3A*VA + y3B*VB + y3C*VC)^2
    r2y3B = -k2*y3B*VB/(y3A*VA + y3B*VB + y3C*VC)^2
    r2y3C = -k2*y3B*VC/(y3A*VA + y3B*VB + y3C*VC)^2

    hout = zeros(typeof(p[1]),11,11)
    hout[1,1] = 1.0
    hout[1,6] = 1.0
    hout[1,2] = -1.0

    hout[2,2] = 1.0
    hout[2,3] = - y3A
    hout[2,7] = - F3 - r1y3A*V
    hout[2,8] = - r1y3B*V
    hout[2,9] = - r1y3C*V

    hout[3,7] = (r1y3A+r2y3A)*V
    hout[3,8] = (r1y3B+r2y3B)*V
    hout[3,9] = (r1y3C+r2y3C)*V

    hout[4,3] = -y3C
    hout[4,7] = r2y3A*V
    hout[4,8] = r2y3B*V
    hout[4,9] = r2y3C*V - y3C

    hout[5,7] = -1.0
    hout[5,8] = -1.0
    hout[5,9] = -1.0

    hout[6,3] = 1.0
    hout[6,4] = -1.0
    hout[6,6] = -1.0

    hout[7,3] = y3A
    hout[7,6] = -F7
    hout[7,7] = F3

    hout[8,3] = y3B
    hout[8,4] = -y4B
    hout[8,10] = F3
    hout[8,11] = -F4

    hout[9,4] = 1.0
    hout[9,5] = -1.0

    hout[10,4] = y4B
    hout[10,10] = F4

    hout[11,4] = y4C
    hout[11,5] = -1.0
    hout[11,11] = F4

    return hout
end

hjrs1(z,p) = ForwardDiff.jacobian(x-> hrs(x[1:11],x[12:14]),vcat(z,p))[1:11,1:11]

function frs(z,p)
    # translate the inputs
    V = p[1]
    k1 = p[2]
    k2 = p[3]
    F1 = z[1]
    F2 = z[2]
    F3 = z[3]
    F4 = z[4]
    F6 = z[5]
    F7 = z[6]
    y3A = z[7]
    y3B = z[8]
    y3C = z[9]
    y4B = z[10]
    y4C = z[11]

    alpha = 2.5
    beta = 0.52

    C_steam = 21.67E-3
    C_cool = 4.65E-3

    C_cap1 = 132718 + F3*(369*y3A - 1113.9*y3B)
    C_cap2 = 25000 + F4*(6984.5*y4B - 3869.53*y4C^2)
    C_capCSTR = 25764 + 8178*V
    C_op1 = F3*(3 + 36.11*(1.0-y4B-y4C) + 7.71*y4B)*(C_steam + C_cool)
    C_op2 = F4*(26.21 + 29.45*y4B)*(C_steam + C_cool)
    C_ann = (1.0/alpha)*(C_cap1 + C_cap2 + C_capCSTR) + beta*(C_op1 + C_op2)
    return C_ann
end

function LBDrs_func(i)
  Yrs = IntervalBox(Interval(50.0,150.0),Interval(100.0,300.0),
                    Interval(100.0,300.0),Interval(50.0,60.0),
                    Interval(0.0,10.0),Interval(40.0,250.0),
                    Interval(0.5,0.9),Interval(0.1,0.5),
                    Interval(0.001,0.2),Interval(0.142857,1.0),
                    Interval(0.142857E-3,1.0),
                    Interval(14.0,24.0),Interval(0.405,0.415),Interval(0.0545,0.0555))
  return Yrs[i].lo
end

function UBDrs_func(i)
       Yrs = IntervalBox(Interval(50.0,150.0),Interval(100.0,300.0),
                         Interval(100.0,300.0),Interval(50.0,60.0),
                         Interval(0.0,10.0),Interval(40.0,250.0),
                         Interval(0.5,0.9),Interval(0.1,0.5),
                         Interval(0.001,0.2),Interval(0.142857,1.0),
                         Interval(0.142857E-3,1.0),
                         Interval(14.0,24.0),Interval(0.405,0.415),Interval(0.0545,0.0555))
       return Yrs[i].hi
end

imp_opt1 = ImplicitSolver()
imp_opt1.opts.style = MC_style #"KrawczykCW" # sets MC contractor style
imp_opt1.h = hrs
imp_opt1.hj = hjrs1
imp_opt1.f = frs
imp_opt1.g = [] # inherites bounds from jump model
imp_opt1.nx = 11
imp_opt1.flag = true
imp_opt1.Intv_Cntr = intv_style
imp_opt1.ParamInt = opt1
imp_opt1.numConstr = 0
imp_opt1.gL_Loc = []
imp_opt1.gU_Loc = []
imp_opt1.gL = []
imp_opt1.gU = []

fxp = frs(ones(11),ones(3))
hxp = hrs(ones(11),ones(3))
hjxp = hjrs(ones(11),ones(3))

jumpmodel4 = Model(solver=EAGO_NLPSolver(ImplicitOpts = imp_opt1,
                                         LBD_func_relax = LBD_fr,
                                         LBD_problem_relax = LBD_pr,
                                         LBD_problem_solver = LBD_ps,
                                         UBD_func_relax = "Original",
                                         UBD_problem_relax = "NLP2",
                                         UBD_problem_solver = "Ipopt",
                                         probe_depth = -1,
                                         variable_depth = 1000,
                                         STD_RR_depth = -1,
                                        DAG_depth = -1))

x = @variable(jumpmodel4, [i=1:14], lowerbound=LBDrs_func(i), upperbound=UBDrs_func(i))

@NLconstraint(jumpmodel4, x[4] + x[9] - x[5] == 0.0 )
@NLconstraint(jumpmodel4, x[5] - x[10]*x[6] - (x[2]*x[10]/(x[10]*(8.937E-2) + x[11]*(1.018E-1) + x[12]*(1.130E-1)))*x[1] == 0.0 )
@NLconstraint(jumpmodel4, ((x[2]*x[10]/(x[10]*(8.937E-2) + x[11]*(1.018E-1) + x[12]*(1.130E-1)))+(x[3]*x[11]/(x[10]*(8.937E-2) +
                            x[11]*(1.018E-1) + x[12]*(1.130E-1))))*x[1] - 50 == 0.0 )

@NLconstraint(jumpmodel4, (x[3]*x[11]/(x[10]*(8.937E-2) + x[11]*(1.018E-1) + x[12]*(1.130E-1)))*x[1] - x[12]*x[6] == 0.0 )
@NLconstraint(jumpmodel4, 1.0 - x[10] - x[11] - x[12] == 0.0 )
@NLconstraint(jumpmodel4, x[6] - x[7] - x[9] == 0.0 )

@NLconstraint(jumpmodel4, x[10]*x[6] - x[9] == 0.0 )
@NLconstraint(jumpmodel4, x[11]*x[6] - x[13]*x[7] == 0.0 )
@NLconstraint(jumpmodel4, x[7] - 50 - x[8] == 0.0 )

@NLconstraint(jumpmodel4, x[13]*x[7] - 50 == 0.0 )
@NLconstraint(jumpmodel4, x[14]*x[7] - x[8] == 0.0 )

@NLobjective(jumpmodel4, Min, (1.0/2.5)*(132718 + x[3]*(369*x[7] - 1113.9*x[8]) + 25000 + x[4]*(6984.5*x[10] - 3869.53*x[11]^2) + 25764 + 8178*x[12]) + 0.52*(x[3]*(3 +
                               36.11*(1.0-x[13]-x[14]) + 7.71*x[11])*(21.67E-3 + 4.65E-3) + x[4]*(26.21 + 29.45*x[10])*(21.67E-3 + 4.65E-3)))
status = solve(jumpmodel4)
=#
