# create the SIP option object for the solver
SIPopt1 = SIPOptions()

# 1D Example (7.4.1 from thesis)
# solution f = -15.8077 @ y = 2.95275
f(x) = (x[1]-3.5)^4 - 5*(x[1]-3.5)^3 - 2*(x[1]-3.5)^2 + 15*(x[1]-3.5)
function h(out,x,y,p)
    out[1] = y[1] - (x[1] - (x[1]^3)/6 + (x[1]^5)/120)/sqrt(y[1]) - p[1]
end
function hj(out,x,y,p)
    out[1] = 1.0 + 0.5*(x[1] - (x[1]^3)/6 + (x[1]^5)/120)/sqrt(y[1])^3
end
gSIP(x,y,p) = y[1] + cos(x[1] - p[1]/90) - p[1]
x_lo, x_hi = [0.5], [8.0]
y_lo, y_hi = [68.8], [149.9]
p_lo, p_hi = [80.0], [120.0]
impout1 = implicit_sip_solve(f, gSIP, h, x_lo, x_hi, y_lo, y_hi,
                             p_lo, p_hi,opts = SIPopt1, hj = hj)

#=
module Check_SIP_Routines

    using Compat
    using Compat.Test
    using EAGO
    #=
    @testset "SIPres parameter errors" begin
        SIPopt1 = SIPOptions()

        function f1(x)
            (1/3)*x[1]^2 + x[2]^2 + x[1]/2
        end
        function gSIP1(x,p)
            (1.0-(x[1]^2)*(p[1]^2))^2 - x[1]*p[1]^2 - x[2]^2 + x[2]
        end
        Xlow = [-1000.0,-1000.]
        Xhigh = [1000.0,1000.0]
        Plow = [0.0]
        Phigh = [1.0]
        SIPoutput1 = explicit_sip_solve(f1, gSIP1, Xlow, Xhigh, Plow, Phigh, opts = SIPopt1)

        SIPopt1 = SIPOptions()

        f1(x) = (1/3)*x[1]^2 + x[2]^2 + x[1]/2
        gSIP1(x,p) = (1.0-(x[1]^2)*(p[1]^2))^2 - x[1]*p[1]^2 - x[2]^2 + x[2]
        Xlow = [-1000.0,-1000.0]
        Xhigh = [1000.0,1000.0]
        Plow = [0.0]
        Phigh = [1.0]

        SIPopt1.r0 = 2.0
        SIPopt1.eps_g0 = -0.9
        @test_throws ErrorException EAGO.explicit_sip_solve(f1, gSIP1, Xlow, Xhigh, Plow, Phigh, opts = SIPopt1)

        SIPopt1.r0 = 0.1
        SIPopt1.eps_g0 = 0.9
        @test_throws ErrorException EAGO.explicit_sip_solve(f1, gSIP1, Xlow, Xhigh, Plow, Phigh, opts = SIPopt1)

        SIPopt1.r0 = 2.0
        SIPopt1.LBP_Opt.DAG_depth = 100
        @test_throws ErrorException EAGO.explicit_sip_solve(f1, gSIP1, Xlow, Xhigh, Plow, Phigh, opts = SIPopt1)
    end
    =#

    @testset "impSIPres parameter errors" begin

              # create the SIP option object for the solver
              SIPopt1 = SIPOptions()

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
              x_lo, x_hi = [0.5], [8.0]
              y_lo, y_hi = [68.8], [149.9]
              p_lo, p_hi = [80.0], [120.0]
              @test_throws ErrorException implicit_sip_solve(f, h, gSIP, x_lo, x_hi, y_lo, y_hi,
                                                             p_lo, p_hi,opts = SIPopt1, hj = hj)

              SIPopt1.r0 = 2.0
              SIPopt1.eps_g0 = -0.9
              @test_throws ErrorException implicit_sip_solve(f, h, gSIP, x_lo, x_hi, y_lo, y_hi,
                                                             p_lo, p_hi,opts = SIPopt1, hj = hj)

              SIPopt1.eps_g0 = 0.9
              SIPopt1.LBP_Opt.DAG_depth = 100
              @test_throws ErrorException implicit_sip_solve(f, h, gSIP, x_lo, x_hi, y_lo, y_hi,
                                                             p_lo, p_hi,opts = SIPopt1, hj = hj)

              SIPopt1.LBP_Opt.DAG_depth = -1
              SIPopt1.UBP_Opt.DAG_depth = 100
              @test_throws ErrorException implicit_sip_solve(f, h, gSIP, x_lo, x_hi, y_lo, y_hi,
                                                             p_lo, p_hi,opts = SIPopt1, hj = hj)
    end

    #=
    @testset "SIPres test problem (SipTestSet #2)" begin
        xl = [-1000.0, -1000.0]
        xu = [1000.0, 1000.0]
        pl = [0.0]
        pu = [1.0]
        f(x) = (x[1]^2)/3.0 + x[2]^2 + x[1]/2.0
        g(x,p) = (1.0 - (x[1]^2)*(p[1]^2)) - x[1]*p[1]^2 - x[2] + x[2]
    end

    @testset "SIPres test problem (SipTestSet #5)" begin
        xl = [-1000.0, -1000.0, -1000.0]
        xu = [1000.0, 1000.0, 1000.0]
        pl = [0.0]
        pu = [1.0]
        f(x) = exp(x[1]) + exp(x[2]) + exp(x[3])
        g(x,p) = 1.0/(1.0 + p[1]^2) - x[1] - x[2]*p[1] - x[3]*p[1]^2
    end

    @testset "SIPres test problem (SipTestSet #4.3)" begin
        xl = [-1000.0, -1000.0, -1000.0]
        xu = [1000.0, 1000.0, 1000.0]
        pl = [0.0]
        pu = [1.0]
        f(x) = x[1] + x[2]/2.0 + x[3]/3.0
        g(x,p) = exp(p[1] - 1.0) - x[1] - x[2]*p[1] - x[3]*p[1]^2
    end

    @testset "SIPres Infeasible problem" begin
        xl = [-1000.0, -1000.0, -1000.0]
        xu = [1000.0, 1000.0, 1000.0]
        pl = [0.0]
        pu = [1.0]
        f(x) = x[1] + x[2]/2.0 + x[3]/3.0
        g(x,p) = x[1] + p[1]^2 + 1100.0
    end
    =#

    @testset "impSIPres test problem (1D Example, 7.4.1 from Stuber thesis)" begin

        # create the SIP option object for the solver
        SIPopt1 = SIPOptions()

        SIPopt1.eps_g0 = 0.9
        SIPopt1.tol = 1E-2              # SIP tolerance
        SIPopt1.r0 = 2.0                # reduction factor for SIP routine
        SIPopt1.kmax = 5                # maximum number of iteration for SIP routine
        SIPopt1.inn_tol = 0.0           # tolerance factor usually set to tolerance of inner program

        # solution f = -15.8077 @ y = 2.95275
        f(x) = (x[1]-3.5)^4 - 5*(x[1]-3.5)^3 - 2*(x[1]-3.5)^2 + 15*(x[1]-3.5)
        function h!(out,x,y,p)
           out[1] = [y[1]-(x[1]-(x[1]^3)/6+(x[1]^5)/120)/sqrt(y[1])-p[1]]
        end
        function hj!(out,x,y,p)
           out[1] = [1.0+(x[1]-(x[1]^3)/6+(x[1]^5)/120)/(2.0*sqrt(y[1]^3))]
        end
        gSIP(x,y,p) = y[1] + cos(x[1]-p[1]/90) - p[1]
        x_lo, x_hi = [0.5], [8.0]
        y_lo, y_hi = [68.8], [149.9]
        p_lo, p_hi = [80.0], [120.0]
        impout1 = implicit_sip_solve(f, h!, gSIP, x_lo, x_hi, y_lo, y_hi,
                                     p_lo, p_hi,opts = SIPopt1, hj = hj!)
        # get solution values
        UBD = impout1.UBD               # upper bound
        LBD = impout1.LBD               # lower bound
        feas = impout1.feas             # is problem feasible?

        @test isapprox(LBD,-7.898552614446456,atol=1E-3)
        @test isapprox(UBD,-7.896382235986764,atol=1E-3)
        @test feas == true
    end

    @testset "impSIPres test problem (impSIPres test problem Stuber2015, Ex2)" begin
        # create the SIP option object for the solver
        SIPopt1 = SIPOptions()

        SIPopt1.eps_g0 = 0.9
        SIPopt1.tol = 1E-2              # SIP tolerance
        SIPopt1.r0 = 2.0                # reduction factor for SIP routine
        SIPopt1.kmax = 5                # maximum number of iteration for SIP routine
        SIPopt1.inn_tol = 0.0           # tolerance factor usually set to tolerance of inner program


        # input data
        A = [7.00961, 7.00877, 6.9895]
        B = [1022.48, 1134.15, 1216.92]
        C = [248.145, 238.678, 227.451]
        z = [0.5, 0.5, 0.1]

        f(x) = -x[2]
        function h!(out,x,y,p)
           K = exp10.(A - B./(C + x[1]))/p[1]
           out[1] = sum((z[i]*K[i]-1.0)/((K[i]-1.0)*y[1]+1.0 for i in 1:3)
        end
        function hj!(out,x,y,p)
           out[1] = # TO ADD
        end
        function gSIP(x,y,p)
            K3 = exp10.(A[3] - B[3]./(C[3] + x[1]))/p[1]
            return x[1] - z[3]*K3/((K3-1.0)*y+1.0)
        end
        x_lo, x_hi = [80.0, -1.0], [90.0, 1.0]
        y_lo, y_hi = [0.0], [0.1]
        p_lo, p_hi = [4400.0], [5100.0]
        impout1 = implicit_sip_solve(f, h!, gSIP, x_lo, x_hi, y_lo, y_hi,
                                     p_lo, p_hi,opts = SIPopt1, hj = hj!)
        # get solution values
        UBD = impout1.UBD               # upper bound
        LBD = impout1.LBD               # lower bound
        feas = impout1.feas             # is problem feasible?

        #@test isapprox(LBD,-$TBD,atol=1E-3)
        #@test isapprox(UBD,$TBD,atol=1E-3)
        @test feas == true
    end

    #=
    @testset "impSIPres test problem (impSIPres test problem Stuber2015, Ex3)" begin
        # create the SIP option object for the solver
        SIPopt1 = SIPOptions()

        SIPopt1.eps_g0 = 0.9
        SIPopt1.tol = 1E-2              # SIP tolerance
        SIPopt1.r0 = 2.0                # reduction factor for SIP routine
        SIPopt1.kmax = 5                # maximum number of iteration for SIP routine
        SIPopt1.inn_tol = 0.0           # tolerance factor usually set to tolerance of inner program


        # input data

        # functional form
        f(x) = x[1]
        function h!(out,x,y,p)
            k1 = p[1]
            k2 = p[2]
            F1 = p[3]
        end
        function hj!(out,x,y,p)
        end
        gSIP(x,y,p) = 22.0 - y[2]*y[4]

        # bounds
        x_lo, x_hi = [10.0], [20.0]
        y_lo, y_hi = [0.15, 0.3, 0.0, 60.0], [0.85, 0.65, 0.12, 70.0]
        p_lo, p_hi = [0.38, 0.053, 60.0], [0.42, 0.058, 70.0]

        impout1 = implicit_sip_solve(f, h!, gSIP, x_lo, x_hi, y_lo, y_hi,
                                     p_lo, p_hi,opts = SIPopt1, hj = hj!)
        # get solution values
        UBD = impout1.UBD               # upper bound
        LBD = impout1.LBD               # lower bound
        feas = impout1.feas             # is problem feasible?

        #@test isapprox(LBD,-$TBD,atol=1E-3)
        #@test isapprox(UBD,$TBD,atol=1E-3)
        @test feas == true
    end
    =#
end
=#
