#workspace()

using IntervalArithmetic
using EAGO

# create the SIP option object for the solver
SIPopt1 = SIP_opts()

# create solver with specified options options for lower level problem
sep1in = EAGO_NLPSolver(LBD_func_relax = "NS-STD-OFF",  # use standard McCormick relaxations
                        LBDsolvertype = "LP",           # use an LP problem structure for relaxed problems
                        UBDsolvertype = "Ipopt",        # use NLP solver upper bounds (currently preferred solver)
                        probe_depth = -1,               # disable probing
                        variable_depth = 1000,          # use duality based range reduction to a depth of 1000 (use to high depth recommended)
                        DAG_depth = -1,                 # don't use a DAG contractor (I need to update this for implicit SIP)
                        STD_RR_depth = -1,              # don't use standard range reduction (problems get quite large)
                        verbosity = "None",             # specify printing level for global optimization problem
                        validated = true,               # use numerically validated intervals
                        atol = 1E-7,                    # absolute tolerance (May need to play with this)
                        rtol = 1E-5)                    # relative tolerance (May need to play with this)

# create a solver for the lower/upper problems
sep1lu = EAGO_NLPSolver(LBD_func_relax = "NS-STD-OFF",
                        LBDsolvertype = "LP",
                        UBDsolvertype = "Ipopt",
                        probe_depth = -1,
                        variable_depth = 1000,
                        DAG_depth = -1,
                        STD_RR_depth = -1,
                        verbosity = "None",
                        validated = true,
                        atol = 1E-7,
                        rtol = 1E-5)

SIPopt1.LLP_Opt = sep1in        # Set solver for use in lower level problem
SIPopt1.LBP_Opt = sep1lu        # Set solver for use in lower bounding problem
SIPopt1.UBP_Opt = sep1lu        # Set solver for use in upper bounding problem

SIPopt1.eps_g0 = 0.9
SIPopt1.tol = 1E-2              # SIP tolerance
SIPopt1.r0 = 2.0                # reduction factor for SIP routine
SIPopt1.kmax = 5                # maximum number of iteration for SIP routine
SIPopt1.inn_tol = 0.0           # tolerance factor usually set to tolerance of inner program

# 1D Example (7.4.1 from thesis)
# solution f = -15.8077 @ y = 2.95275
f(x) = (x[1]-3.5)^4 - 5*(x[1]-3.5)^3 - 2*(x[1]-3.5)^2 + 15*(x[1]-3.5)
function h(x,y,p)
    [y[1]-(x[1]-(x[1]^3)/6+(x[1]^5)/120)/sqrt(y[1])-p[1]]
end
function hj(x,y,p)
    [1.0+(x[1]-(x[1]^3)/6+(x[1]^5)/120)/(2.0*sqrt(y[1]^3))]
end
gSIP(x,y,p) = y[1] + cos(x[1]-p[1]/90) - p[1]
xBnds = [Interval(0.5,8.0)]
yBnds = [Interval(68.8,149.9)]
pBnds = [Interval(80,120)]
impout1 = Implicit_SIP_Solve(f,h,hj,gSIP,xBnds,yBnds,pBnds,SIPopt1)

# get solution values
k = impout1.k                   # number of iterations
UBD = impout1.UBD               # upper bound
LBD = impout1.LBD               # lower bound
feas = impout1.feas             # is problem feasible?
LBP_time = impout1.LBP_time     # time spent solving lower bounding problem
LLP_time = impout1.LLP_time     # time spent solving lower level problem
UBP_time = impout1.UBP_time     # time spent solving upper bounding problem
xbar = impout1.xbar             # solution value
