module ImplicitSIP_Tests

using Compat
using Compat.Test
using IntervalArithmetic
using EAGO

@testset "Test Reformulation Functions" begin

gSIP(x,y,p) = [x[1]^2+y[1]-p[2]]
g(x,y,p) = [x[1]+y[1]+p[1], 2.0x[1]+3.0y[2]+4.0p[2]]
h(x,y,p) = [x[1]+x[2]+y[1]-p[2];
            x[1]-x[2]-y[2]+p[1]]
hj(x,y,p) = [1.0 1.0;
             1.0 -y[2]]
Pset = [[1.0; 2.0],[3.0; 4.0]]
eps_g = 1.5
ny = 2
ng = 2

X = [Interval(1,2);Interval(3,4)]
Y = [Interval(4,5);Interval(5,6)]
x = [4.5,5.5]
y = [1.5,3.5]
Ylo1,Yhi1,st1,p1 = EAGO.Reform_Imp_Y(X,Y,Pset)

@test Ylo1 == [4.0, 5.0, 4.0, 5.0, 1.0, 3.0]
@test Yhi1 == [5.0, 6.0, 5.0, 6.0, 2.0, 4.0]
@test st1 == 4
@test p1 == 6

yref = [1.5,3.5,1.5,3.5]
out3 = EAGO.Reform_Imp_H(h,x,yref,Pset,ny)
@test out3 == [9.5, -3.5, 7.5, -1.5]

out4 = EAGO.Reform_Imp_G(g,x,yref,Pset,ny,eps_g)
@test out4 == [8.5, 29.0, 10.5, 37.0]

yrefa = [1.5,3.5,1.5,2.5]
out6 = EAGO.Reform_Imp_HJ(hj,x,yref,Pset,ny)
out6a = EAGO.Reform_Imp_HJ(hj,x,yrefa,Pset,ny)
@test out6 ==  [1.0   1.0  0.0   0.0
          1.0  -3.5  0.0   0.0
          0.0   0.0  1.0   1.0
          0.0   0.0  1.0  -3.5]
@test out6a ==  [1.0   1.0  0.0   0.0
          1.0  -3.5  0.0   0.0
          0.0   0.0  1.0   1.0
          0.0   0.0  1.0  -2.5]

# later if necessary
out5 = EAGO.Reform_Imp_HG(h,g,x,yrefa,Pset,ny,ng,eps_g)
@test out5 == [11.0, -2.0, 9.0, 1.0, -11.0, 2.0, -9.0, -1.0, 8.5, 29.0, 10.5, 34.0]

end

#=
# solves example SIP #1 with DAG contractor disabled
@testset "SemiInfinite Explicit RHS" begin
SIPopt1 = SIP_opts()
sep1lu = EAGO_NLPSolver(LBD_func_relax = "NS-STD-OFF",
                        LBDsolvertype = "LP",
                        probe_depth = -1,
                        variable_depth = 1000,
                        DAG_depth = -1,
                        STD_RR_depth = 1000,
                        ImplicitFlag = true,
                        verbosity = "Normal",
                        validated = true))
sep1lu.BnBSolver.Verbosity = "None"
sep1in = EAGO_NLPSolver(BD_func_relax = "NS-STD-OFF",
                            LBDsolvertype = "LP",
                            probe_depth = -1,
                            variable_depth = 1000,
                            DAG_depth = -1,
                            STD_RR_depth = 1000,
                            ImplicitFlag = true,
                            verbosity = "Normal",
                            validated = true))

sep1in.BnBSolver.Verbosity = "None"
SIPopt1.LLP_Opt = sep1in
SIPopt1.LBP_Opt = sep1lu
SIPopt1.UBP_Opt = sep1lu
f1(x) = (1/3)*x[1]^2 + x[2]^2 + x[1]/2
gSIP1(x,p) = (1.0-(x[1]^2)*(p[1]^2))^2 - x[1]*p[1]^2 - x[2]^2 + x[2]
X1 = [MCInterval(-1000.0,1000.0),MCInterval(-1000.0,1000.0)]
P1 = [MCInterval(0.0,1.0)]
SIPoutput1 = Implicit_SIP_Solve(f1,gSIP1,X1,P1,SIPopt1)

@test isapprox(SIPoutput1.LBD,0.19452787006676814,atol=1E-3)
@test isapprox(SIPoutput1.UBD,0.19452787006676814,atol=1E-3)

end
=#
end
