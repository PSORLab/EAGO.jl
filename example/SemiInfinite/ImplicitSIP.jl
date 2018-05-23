#workspace()

using IntervalArithmetic
using EAGO

# create the SIP option object for the solver
SIPopt1 = SIP_opts()

# create solver with specified options options for lower level problem
sep1in = EAGO_NLPSolver(LBD_func_relax = "NS-STD-OFF",  # use standard McCormick relaxations
                        LBDsolvertype = "LP",           # use an LP problem structure for relaxed problems
                        #UBDsolvertype = "Interval",
                        UBDsolvertype = "Ipopt",
                        probe_depth = -1,               # disable probing
                        variable_depth = -1000,          # use duality based range reduction to a depth of 1000
                        DAG_depth = -1,                 # don't use a DAG contractor
                        STD_RR_depth = -1,            # use standard range reduciton to a depth of 1000
                        verbosity = "None",           # specify printing level for global optimization problem
                        validated = true,             # use numerically validated intervals
                        atol = 1E-7,
                        rtol = 1E-5)

# create a solver for the lower/upper problems
sep1lu = EAGO_NLPSolver(LBD_func_relax = "NS-STD-OFF",
                        LBDsolvertype = "LP",
                        UBDsolvertype = "Ipopt",
                        probe_depth = -1,
                        variable_depth = -1000,
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
SIPopt1.tol = 1E-4
SIPopt1.r0 = 2.0
SIPopt1.kmax = 4

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

#=
# Flash Example
# Solution: f = y = 10.1794m^3  @ p = (0.38 0.058 60) T
# defines implicit function
function hCSTR!(hout,z,y,p)

    VA = 8.397E-2
    VB = 1.018E-1
    VC = 1.13E-1

    r1 = p[1]*z[1]/(z[1]*VA+z[2]*VB+z[3]*VC)
    r2 = p[2]*z[2]/(z[1]*VA+z[2]*VB+z[3]*VC)

    zA1 = 1.0
    zA2 = 0.0
    zA3 = 0.0

    hout[1] = p[3] - z[1]*z[4] - y[1]*r[1]
    hout[2] = - z[2]*z[4] - y[1]*(r[1]-r[2])
    hout[3] = - z[3]*z[4] - y[1]*r[2]
    hout[4] = 1 - z[1] - z[2] - z[3]
end
function hjCSTR!(hout,z,y,p)

        VA = 8.397E-2
        VB = 1.018E-1
        VC = 1.13E-1

        r1 = p[1]*z[1]/(z[1]*VA+z[2]*VB+z[3]*VC)
        r1P1 = z[1]/(z[1]*VA+z[2]*VB+z[3]*VC)
        r1Z1 = p[1]*(VC*z[3]+VB*z[2])/(z[1]*VA+z[2]*VB+z[3]*VC)^2
        r1Z2 = -p[1]*VB*z[1]/(z[1]*VA+z[2]*VB+z[3]*VC)^2
        r1Z3 = -p[1]*VC*z[1]/(z[1]*VA+z[2]*VB+z[3]*VC)^2

        r2 = p[2]*z[2]/(z[1]*VA+z[2]*VB+z[3]*VC)
        r2P1 = z[2]/(z[1]*VA+z[2]*VB+z[3]*VC)
        r2Z1 = -p[2]*VB*z[2]/(z[1]*VA+z[2]*VB+z[3]*VC)^2
        r2Z2 = p[2]*(VC*z[3]+VA*z[1])/(z[1]*VA+z[2]*VB+z[3]*VC)^2
        r2Z3 = -p[2]*VC*z[2]/(z[1]*VA+z[2]*VB+z[3]*VC)^2

        zA1 = 1.0
        zA2 = 0.0
        zA3 = 0.0

        hout[1,1] = -z[4] - y[1]*r1Z1 - r1
        hout[1,2] = -y[1]*r1Z2
        hout[1,3] = -y[1]*r1Z3
        hout[1,4] = -z[1]
        hout[1,5] = -y[1]*r1P1
        hout[1,7] = 1.0

        hout[2,1] = - y[1]*(r1Z1-r2Z1) - (r[1]-r[2])
        hout[2,2] = - z[4] - y[1]*(r1Z2-r2Z2)
        hout[2,3] = - y[1]*(r1Z3-r2Z3)
        hout[2,4] = - z[2]
        hout[2,5] = - y[1]*(r1P1)
        hout[2,6] = - y[1]*(-r2P1)

        hout[3,1] = - r[2] - y[1]*r2Z1
        hout[3,2] = - y[1]*r2Z2
        hout[3,3] = - z[4] - y[1]*r2Z3
        hout[3,4] = - z[3]
        hout[3,6] = - y[1]*r[2]

        hout[4,1] = -1.0
        hout[4,2] = -1.0
        hout[4,3] = -1.0
end
fCSTR(z,y,p) = y[1]
gSIPCSTR(z,y,p) = 22-z[2]*z[4]
# defines
YCSTR = [10..20]
PCSTR = [0.38..0.42,0.053..0.058,60..70]
XCSTR = [0.15..0.85,0.3..0.65,0.0..0.12,60..70]
impout2 = Implicit_SIP_Solve(fCSTR,hCSTR!,hjCSTR!,gSIPCSTR,XCSTR,YCSTR,PCSTR,SIPopt1,4)
=#
