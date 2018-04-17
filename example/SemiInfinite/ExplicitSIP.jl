workspace()

using IntervalArithmetic
using EAGO

# solves example SIP #1 with DAG contractor disabled
SIPopt1 = SIP_opts()
sep1lu = EAGO_NLPSolver(probe_depth = -1,
                        variable_depth = 1000,
                        DAG_depth = -1,
                        STD_RR_depth = -1,
                        UBDsolvertype= "Ipopt")
sep1lu.BnBSolver.Verbosity = "None"
sep1in = EAGO_NLPSolver(probe_depth = -1,
                        variable_depth = 1000,
                        DAG_depth = -1,
                        STD_RR_depth = -1,
                        UBDsolvertype= "Ipopt")
sep1in.BnBSolver.Verbosity = "None"
SIPopt1.LLP_Opt = sep1in
SIPopt1.LBP_Opt = sep1lu
SIPopt1.UBP_Opt = sep1lu
f1(x) = (1/3)*x[1]^2 + x[2]^2 + x[1]/2
gSIP1(x,p) = (1.0-(x[1]^2)*(p[1]^2))^2 - x[1]*p[1]^2 - x[2]^2 + x[2]
X1 = [MCInterval(-1000.0,1000.0),MCInterval(-1000.0,1000.0)]
P1 = [MCInterval(0.0,1.0)]
SIPoutput1 = Explicit_SIP_Solve(f1,gSIP1,X1,P1,SIPopt1)

#=
# solves example SIP #1 with DAG contractor enabled
SIPopt2 = SIP_opts()
sep2lu = EAGO_NLPSolver(probe_depth = -1,
                        variable_depth = 1000,
                        DAG_depth = -1,
                        STD_RR_depth = -1)
sep2lu.BnBSolver.Verbosity = "None"
sep2in = EAGO_NLPSolver(probe_depth = -1,
                        variable_depth = 1000,
                        DAG_depth = -1,
                        STD_RR_depth = -1)
sep2in.BnBSolver.Verbosity = "None"
SIPopt2.LLP_Opt = sep1in
SIPopt2.LBP_Opt = sep1lu
SIPopt2.UBP_Opt = sep1lu
f2(x) = (1/3)*x[1]^2 + x[2]^2 + x[1]/2
gSIP2(x,p) = (1-(x[1]^2)*(p[1]^2))^2 - x[1]*p[1]^2 - x[2]^2 + x[2]
# note that p variable is moved to end of x variable array
gSIP2ex = :((1.0-(x[1]^2)*(x[3]^2))^2 - x[1]*x[3]^2 - x[2]^2 + x[2])
SIPopt2.gSIPExp = gSIP2ex
X2 = [MCInterval(-1000.0,1000.0),MCInterval(-1000.0,1000.0)]
P2 = [MCInterval(0.0,1.0)]
SIPoutput2 = Explicit_SIP_Solve(f2,gSIP2,X2,P2,SIPopt2)

# solves example SIP #7 with data contractor disabled
SIPopt3 = SIP_opts() # SIP parameter options for file
sep3lu = EAGO_NLPSolver()  # solver selection for lower/upper bounding problems
sep3lu.BnBSolver.Verbosity = "None" # lower/upper bounding problems console output level
sep3in = EAGO_NLPSolver() # solver selection for inner problem
sep3in.BnBSolver.Verbosity = "None" # inner problem console output level
SIPopt3.LLP_Opt = sep3in # set inner problem solver
SIPopt3.LBP_Opt = sep3lu # set lower problem solver
SIPopt3.UBP_Opt = sep3lu # set upper problem solver

f3(x) = x[1]^2 + x[2]^2 + x[3]^2
gSIP3(x,p) = x[1]*(p[1]+p[2]^2+1.0)+x[2]*(p[1]*p[2]-p[2]^2)+x[3]*(p[1]*p[2]+p[2]^2+p[2])+1.0
X3 = [MCInterval(-1000,1000),MCInterval(-1000,1000),MCInterval(-1000,1000)]
P3 = [MCInterval(0,1),MCInterval(0,1)]
SIPoutput2 = Explicit_SIP_Solve(f3,gSIP3,X3,P3,SIPopt3)
=#
