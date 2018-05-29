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
          @test out5 == [11.0, -2.0, 9.0, 1.0, 8.5, 29.0, 10.5, 34.0]

end

# solves example SIP #1 with DAG contractor disabled
@testset "SemiInfinite Implicit RHS" begin

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
          UBD = impout1.UBD               # upper bound
          LBD = impout1.LBD               # lower bound
          feas = impout1.feas             # is problem feasible?

          @test isapprox(LBD,-7.898552614446456,atol=1E-3)
          @test isapprox(UBD,-7.896382235986764,atol=1E-3)
          @test feas == true
end

@testset "SemiInfinite Error Handling" begin

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
          SIPopt1.r0 = 0.1                # reduction factor for SIP routine
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
          @test_throws ErrorException Implicit_SIP_Solve(f,h,hj,gSIP,xBnds,yBnds,pBnds,SIPopt1)

          SIPopt1.r0 = 2.0
          SIPopt1.eps_g0 = -0.9
          @test_throws ErrorException Implicit_SIP_Solve(f,h,hj,gSIP,xBnds,yBnds,pBnds,SIPopt1)

          SIPopt1.eps_g0 = 0.9
          SIPopt1.LBP_Opt.DAG_depth = 100
          @test_throws ErrorException Implicit_SIP_Solve(f,h,hj,gSIP,xBnds,yBnds,pBnds,SIPopt1)

          SIPopt1.LBP_Opt.DAG_depth = -1
          SIPopt1.UBP_Opt.DAG_depth = 100
          @test_throws ErrorException Implicit_SIP_Solve(f,h,hj,gSIP,xBnds,yBnds,pBnds,SIPopt1)
end
