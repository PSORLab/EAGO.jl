# package setup

#workspace()
using EAGO
using Ipopt
using JuMP
using StaticArrays
using IntervalArithmetic

########################## JuMP+ Interface #####################################

# Solves the quadratically constrained problem (Example 5.1, Stuber2015)
h1(x,p) = [x[1]^2 + p[1]*x[1] + 4.0]
hj1(x,p) = [2.0*x[1] + p[1]]
f1(x,p) = x[1]
LBD1a_func(i) = (i==1) ? (-0.78) : (6.0)
UBD1a_func(i) = (i==1) ? (-0.4) : (9.0)
LBD1b_func(i) = (i==1) ? (-10.0) : (6.0)
UBD1b_func(i) = (i==1) ? (-5.0) : (9.0)
jm1a = Model(solver=EAGO_NLPSolver(LBD_func_relax = "NS-STD",
                                   LBDsolvertype = "LP",
                                   variable_depth = 1000,
                                   ImplicitFlag = true,
                                   verbosity = "Normal",
                                   validated = true))
xa = @variable(jm1a, [i=1:2], lowerbound=LBD1a_func(i), upperbound=UBD1a_func(i))
@NLconstraint(jm1a, xa[1]^2 + xa[2]*xa[1] + 4 == 0.0 )
@NLobjective(jm1a, Min, xa[1])
status1a = Solve_Implicit(jm1a,f1,h1,hj1,x->[],1)


jm1b = Model(solver=EAGO_NLPSolver(LBD_func_relax = "NS-STD",
                                   LBDsolvertype = "LP",
                                   variable_depth = 1000,
                                   ImplicitFlag = true,
                                   verbosity = "Normal",
                                   validated = true))
xb = @variable(jm1b, [i=1:2], lowerbound=LBD1b_func(i), upperbound=UBD1b_func(i))
@NLconstraint(jm1b, xb[1]^2 + xb[2]*xb[1] + 4 == 0.0 )
@NLobjective(jm1b, Min, xb[1])
status1b = Solve_Implicit(jm1b,f1,h1,hj1,x->[],1)

########################## NonJuMP Interface #####################################
nx = 1                                 # number of state variables
h1(x,p) = [x[1]^2 + p[1]*x[1] + 4.0]   # equality constaint yeild implicit variable
hj1(x,p) = [2.0*x[1] + p[1]]           # jacobian of prior line
f1(x,p) = x[1]                         # objective function in f(x,p)
fe1(y) = f1(y[1:nx],y[(nx+1):end])     # f(y) where y from 1:nx is x, and (nx+1):(nx+np) is p
he1(y) = h1(y[1:nx],y[(nx+1):end])     # f(y) where y from 1:nx is x, and (nx+1):(nx+np) is p

# Build the explicit problem model
m = Build_Script(fe1,           # explicitly defined objective
                 [-0.78,6.0],   # lower bounds on y
                 [-0.4,9.0],    # upper bounds on y
                 h = he1,
                 solver = EAGO_NLPSolver(LBD_func_relax = "NS-STD",
                                         LBDsolvertype = "LP",
                                         UBDsolvertype = "Ipopt",
                                         variable_depth = 1000,
                                         ImplicitFlag = true,
                                         verbosity = "Normal",
                                         validated = true))

# Adds the implicitly defined expressions
EAGO.set_Implicit_Model!(m,
                    f1,    # implicitly defined objective
                    h1,    # implicitly defined equality constrain
                    hj1,
                    x->[],
                    nx)

MathProgBase.optimize!(m)            # Solves the model
status = MathProgBase.status(m)      # Returns feasibility status
obj = MathProgBase.getobjval(m)      # Returns objective value
sol = MathProgBase.getsolution(m)    # Returns solution
