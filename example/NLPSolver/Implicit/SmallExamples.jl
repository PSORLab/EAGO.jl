# package setup

workspace()
using EAGO
using Ipopt
using JuMP

# Solves the quadratically constrained problem (Example 5.1, Stuber2015)
h1(x,p) = [x[1]^2 + p[1]*x[1] + 4.0]
hj1(x,p) = [2.0*x[1] + p[1]]
f1(x,p) = x[1]
LBD1a_func(i) = (i==1) ? (-0.78) : (6.0)
UBD1a_func(i) = (i==1) ? (-0.4) : (9.0)
LBD1b_func(i) = (i==1) ? (-10.0) : (6.0)
UBD1b_func(i) = (i==1) ? (-5.0) : (9.0)

jm1a = Model(solver=EAGO_NLPSolver(LBD_func_relax = "NS-STD-OFF",
                                   LBDsolvertype = "LP",
                                   probe_depth = -1,
                                   variable_depth = 1000,
                                   DAG_depth = -1,
                                   STD_RR_depth = 1000,
                                   ImplicitFlag = true,
                                   verbosity = "Normal",
                                   validated = false))
xa = @variable(jm1a, [i=1:2], lowerbound=LBD1a_func(i), upperbound=UBD1a_func(i))
@NLconstraint(jm1a, xa[1]^2 + xa[2]*xa[1] + 4 == 0.0 )
@NLobjective(jm1a, Min, xa[1])
status1a = Solve_Implicit(jm1a,f1,h1,hj1,x->[],1)


jm1b = Model(solver=EAGO_NLPSolver(LBD_func_relax = "NS-STD-OFF",
                                   LBDsolvertype = "LP",
                                   probe_depth = -1,
                                   variable_depth = 1000,
                                   DAG_depth = -1,
                                   STD_RR_depth = 1000,
                                   ImplicitFlag = true,
                                   verbosity = "Normal",
                                   validated = false))
xb = @variable(jm1b, [i=1:2], lowerbound=LBD1b_func(i), upperbound=UBD1b_func(i))
@NLconstraint(jm1b, xb[1]^2 + xb[2]*xb[1] + 4 == 0.0 )
@NLobjective(jm1b, Min, xb[1])
status1b = Solve_Implicit(jm1b,f1,h1,hj1,x->[],1)

LBD1c_func(i) = (i==1) ? (68.8) : (0.5)
UBD1c_func(i) = (i==1) ? (149.9) : (8.0)
g2(x,p) = x[1] + cos(p[1]-80/90) - 80
function h2(y,x)
    [y[1]-(x[1]-(x[1]^3)/6+(x[1]^5)/120)/sqrt(y[1])-80]
end
function hj2(y,x)
    [1.0+(x[1]-(x[1]^3)/6+(x[1]^5)/120)/(2.0*sqrt(y[1]^3))]
end
f2(x,p) = (p[1]-3.5)^4 - 5*(p[1]-3.5)^3 - 2*(p[1]-3.5)^2 + 15*(p[1]-3.5)


jm1ca = Model(solver=EAGO_NLPSolver(LBD_func_relax = "NS-STD-OFF",
                                   LBDsolvertype = "LP",
                                   UBDsolvertype = "Ipopt",
                                   #LBD_func_relax = "Interval",
                                   #LBDsolvertype = "Interval",
                                   probe_depth = -1,
                                   variable_depth = -1000,
                                   DAG_depth = -1,
                                   STD_RR_depth = -1000,
                                   ImplicitFlag = false,
                                   verbosity = "Normal",
                                   validated = true))
@variable(jm1ca, 68.8 <= a <= 149.9)
@variable(jm1ca, 0.5 <= b <= 8.0)
println("ran me 1")
@NLconstraint(jm1ca, a + cos(b-80/90) - 80 <= 0.0 )
println("ran me 2")
@NLconstraint(jm1ca, a-(b-(b^3)/6+(b^5)/120)/sqrt(a)-80 == 0.0 )
println("ran me 3")
@NLobjective(jm1ca, Min, (b-3.5)^4 - 5*(b-3.5)^3 - 2*(b-3.5)^2 + 15*(b-3.5))
println("ran me 4")
status1b = solve(jm1ca)
#=
jm1c = Model(solver=EAGO_NLPSolver(LBD_func_relax = "NS-STD-OFF",
                                   LBDsolvertype = "LP",
                                   UBDsolvertype = "Ipopt",
                                   #LBD_func_relax = "Interval",
                                   #LBDsolvertype = "Interval",
                                   probe_depth = -1,
                                   variable_depth = -1000,
                                   DAG_depth = -1,
                                   STD_RR_depth = -1000,
                                   ImplicitFlag = false,
                                   verbosity = "Normal",
                                   validated = true))
xc = @variable(jm1c, [i=1:2], lowerbound=LBD1c_func(i), upperbound=UBD1c_func(i))
println("ran me 1")
@NLconstraint(jm1c, xc[1] + cos(xc[2]-80/90) - 80 <= 0.0 )
println("ran me 2")
@NLconstraint(jm1c, 0.0 <= xc[1]-(xc[2]-(xc[2]^3)/6+(xc[2]^5)/120)/sqrt(xc[1])-80 <= 0.0 )
println("ran me 3")
@NLobjective(jm1c, Min, (xc[2]-3.5)^4 - 5*(xc[2]-3.5)^3 - 2*(xc[2]-3.5)^2 + 15*(xc[2]-3.5))
println("ran me 4")
status1b = solve(jm1c)
#status1b = Solve_Implicit(jm1c,f2,h2,hj2,g2,1)
=#

#=
# Solves Kolev-based problem (Example 5.2, Stuber 2015)
function h2(x,p)
    [(1.00*10.0^(-9))*(exp(38x[1])-1)+p[1]*x[1]-1.6722*x[2]+0.6689*x[3]-8.0267
     (1.98*10.0^(-9))*(exp(38x[2])-1)+0.6622*x[1]+p[2]*x[2]+0.6622*x[3]+4.0535
     (1.00*10.0^(-9))*(exp(38x[3])-1)+x[1]-x[2]+p[3]*x[3]-6.0]
end
function hj2(x,p)
    [(38.00*10.0^(-9))*(exp(38x[1]))+p[1]    -1.6722                                 0.6689;
      0.6622                                (38*1.98*10.0^(-9))*(exp(38x[2]))+p[2]         0.6622;
      one(p[1])                           (-one(p[1]))                                 (38.00*10.0^(-9))*exp(38x[3])+p[3]]
end
function f2(x1,x2,x3,p1,p2,p3)
    a = [37.3692 18.5805 6.25]
    c = [0.602 1.211 3.6]
    x = [x1;x2;x3]
    p = [p1,p2,p3]
    temp1 = zero(p[1])
    for j = 1:3
        temp2 = a[j]*(p[j]-c[j])^2
        for i = 1:3
            if i != j
                temp2 += a[i]*(p[i]-c[i])-5.0*((j-1)*(j-2)*(x[2]-x[1])+x[1]-x[2]+x[3])
            end
        end
        temp1 += temp2^2
    end
    return temp1
end
function f2p(x,p)
    a = [37.3692 18.5805 6.25]
    c = [0.602 1.211 3.6]
    temp1 = zero(p[1])
    for j = 1:3
        temp2 = a[j]*(p[j]-c[j])^2
        for i = 1:3
            if i != j
                temp2 += a[i]*(p[i]-c[i])-5.0*((j-1)*(j-2)*(x[2]-x[1])+x[1]-x[2]+x[3])
            end
        end
        temp1 += temp2^2
    end
    return temp1
end
function LBD2a_func(i)
    (i==1) && (return 0.5180)
    (i==2) && (return -3.978)
    (i==3) && (return 0.3296)
    (i==4) && (return 0.703)
    (i==5) && (return 1.436)
    (i==6) && (return 3.61)
end
function UBD2a_func(i)
    (i==1) && (return 0.5847)
    (i==2) && (return -3.0464)
    (i==3) && (return 0.5827)
    (i==4) && (return 0.704)
    (i==5) && (return 1.437)
    (i==6) && (return 3.62)
end
contractor = :Newton
mc_opts1 = mc_opts{Float64}(0.5*one(Float64),2,:Dense,contractor,0,0,(1E-12)*one(Float64))
PIntvParams1 = PIntvParams(:Dense,contractor,1E-6,1E-6,3,1,100)
jm2 = Model(solver=EAGO_NLPSolver(LBD_func_relax = "NS-STD-OFF",
                                  LBDsolvertype = "LP",
                                  probe_depth = -1,
                                  variable_depth = 1000,
                                  DAG_depth = -1,
                                  STD_RR_depth = -1,
                                  ImplicitFlag = true,
                                  verbosity = "Normal",
                                  validated = true,
                                  PSmcOpt = mc_opts1,
                                  PIntOpt = PIntvParams1))
xc = @variable(jm2, [i=1:6], lowerbound=LBD2a_func(i), upperbound=UBD2a_func(i))
@NLconstraint(jm2, (1.00*10^(-9))*(exp(38xc[1])-1)+xc[4]*xc[1]-1.6722*xc[2]+0.6689*xc[3]-8.0267 == 0.0 )
@NLconstraint(jm2, (1.98*10^(-9))*(exp(38xc[2])-1)+0.6622*xc[1]+xc[5]*xc[2]+0.6622*xc[3]+4.0535 + 4.0 == 0.0 )
@NLconstraint(jm2, (1.00*10^(-9))*(exp(38xc[3])-1)+xc[1]-xc[2]+xc[6]*xc[3]-6.0 == 0.0 )

#registerEAGO(jm2, :f2, 6, f2, ad=true)
a = [37.3692 18.5805 6.25]
c = [0.602 1.211 3.6]
@NLobjective(jm2, Min, ((a[1]*(xc[4]-c[1]))^2 + a[2]*(xc[5]-c[2]) + a[3]*(xc[6]-c[3]) - 10.0*(xc[1]-xc[2]+xc[3]) +
                       (a[2]*(xc[5]-c[2]))^2 + a[1]*(xc[4]-c[1]) + a[3]*(xc[6]-c[3]) - 10.0*(xc[1]-xc[2]+xc[3]) +
                       (a[3]*(xc[6]-c[3]))^2 + a[1]*(xc[4]-c[1])-10.0*(-xc[1]+xc[2]+xc[3]) + a[3]*(xc[6]-c[3]))^2)
status1b = Solve_Implicit(jm2,f2p,h2,hj2,x->[],3)
# Target = 626.565
=#

# TEST SNOPT IMPLICIT SOLVER


# TEST IPOPT IMPLICIT SOLVER
