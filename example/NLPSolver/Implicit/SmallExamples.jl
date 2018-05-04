# package setup

#workspace()
using EAGO
using Ipopt
using JuMP

# TESTS LP IMPLICIT SOLVER

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


# Solves Kolev-based problem (Example 5.2, Stuber 2015)
function h2(x,p)
    [(1.00*10^(-9))*(exp(38x[1])-1)+p[1]*x[1]-1.6722*x[2]+0.6689*x[3]-8.0267
     (1.98*10^(-9))*(exp(38x[2])-1)+0.6622*x[1]+p[2]*x[2]+0.6622*x[3]+4.0535
     (1.00*10^(-9))*(exp(38x[3])-1)+x[1]-x[2]+p[3]*x[3]-6.0]
end
function hj2(x,p)
    [(38.00*10^(-9))*(exp(38x[1]))+p[1]    p[1]-1.6722                                 0.6689;
      p[2]                                (38*1.98*10^(-9))*(exp(38x[2]))+p[2]         0.6622;
      one(p[1])                           (-one(p[1]))                                 (38.00*10^(-9))*exp(38x[3])+p[3]]
end
function f2(x,p)
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
    (i==1) && (return -5.0)
    (i==2) && (return -5.0)
    (i==3) && (return -5.0)
    (i==4) && (return 0.6020)
    (i==5) && (return 1.2110)
    (i==6) && (return 3.6)
end
function UBD2a_func(i)
    (i==1) && (return 5.0)
    (i==2) && (return 5.0)
    (i==3) && (return 5.0)
    (i==4) && (return 0.7358)
    (i==5) && (return 1.4801)
    (i==6) && (return 4.4)
end

jm2 = Model(solver=EAGO_NLPSolver(LBD_func_relax = "NS-STD-OFF",
                                  LBDsolvertype = "LP",
                                  probe_depth = -1,
                                  variable_depth = 1000,
                                  DAG_depth = -1,
                                  STD_RR_depth = 1000,
                                  ImplicitFlag = true,
                                  verbosity = "Normal",
                                  validated = false))
xc = @variable(jm2, [i=1:6], lowerbound=LBD1b_func(i), upperbound=UBD1b_func(i))
@NLconstraint(jm2, (1.00*10^(-9))*(exp(38xc[1])-1)+xc[4]*xc[1]-1.6722*xc[2]+0.6689*xc[3]-8.0267 == 0.0 )
@NLconstraint(jm2, (1.98*10^(-9))*(exp(38xc[2])-1)+0.6622*xc[1]+xc[5]*xc[2]+0.6622*xc[3]+4.0535 + 4 == 0.0 )
@NLconstraint(jm2, (1.00*10^(-9))*(exp(38xc[3])-1)+xc[1]-xc[2]+xc[6]*xc[3]-6.0 == 0.0 )
@NLobjective(jm1b, Min, xc[1])
status1b = Solve_Implicit(jm2,f2,h2,hj2,x->[],3)

# TEST SNOPT IMPLICIT SOLVER


# TEST IPOPT IMPLICIT SOLVER
