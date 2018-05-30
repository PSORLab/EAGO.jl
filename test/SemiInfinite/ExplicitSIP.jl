module ExplicitSIP_Tests

using Compat
using Compat.Test
using IntervalArithmetic
using EAGO

@testset "SemiInfinite Utilities" begin
X = EAGO.EAGO_NLPSolver()
sip1 = SIP_opts(X)
@test sip1.tol == 1E-3
end

# solves example SIP #1 with DAG contractor disabled
@testset "SemiInfinite Explicit RHS" begin
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

                        @test isapprox(SIPoutput1.LBD,0.19452787006676814,atol=1E-3)
                        @test isapprox(SIPoutput1.UBD,0.19452787006676814,atol=1E-3)

end

@testset "SemiInfinite Error Handling RHS" begin
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

                        println("error test 1")
                        SIPopt1.r0 = 2.0
                        SIPopt1.eps_g0 = -0.9
                        @test_throws ErrorException Explicit_SIP_Solve(f1,gSIP1,X1,P1,SIPopt1)

                        println("error test 2")
                        SIPopt1.r0 = 0.1
                        SIPopt1.eps_g0 = 0.9
                        @test_throws ErrorException Explicit_SIP_Solve(f1,gSIP1,X1,P1,SIPopt1)

                        println("error test 3")
                        SIPopt1.r0 = 2.0
                        SIPopt1.LBP_Opt.DAG_depth = 100
                        @test_throws ErrorException Explicit_SIP_Solve(f1,gSIP1,X1,P1,SIPopt1)

end

end
