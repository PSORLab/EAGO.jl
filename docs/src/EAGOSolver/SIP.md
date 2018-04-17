## Solving Semi-Infinite Programs

```@docs
Explicit_SIP_Solve
```

## Example of solving a SIP without equality constraints
```julia
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
```

## Solving Semi-Infinite Programs with Equality Constraints
