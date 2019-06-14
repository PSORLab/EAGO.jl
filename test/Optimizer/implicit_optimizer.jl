module Check_Implicit_Optimizer

    using Compat
    using Compat.Test
    using EAGO, MathOptInterface, CPLEX
    const MOI = MathOptInterface

    @testset "1D Implicit Lower Evaluator" begin
        low_eval =  EAGO.ImplicitLowerEvaluator{1}()
        new_node = EAGO.NodeBB(Float64[-0.78,6.0],Float64[-0.40,9.0],-Inf,Inf,0,-1,false)
        ytest = [-0.59, 7.5]; ytest2 = [-0.59,6.5]
        EAGO.set_current_node!(low_eval, new_node)
        @test low_eval.current_node.lower_variable_bounds == Float64[-0.78,6.0]
        @test low_eval.current_node.upper_variable_bounds == Float64[-0.40,9.0]

        f(x,p) = x[1]
        g(x,p) = [x[1]]
        function h(out,x,p)
            t1 = x[1]^2
            t2 = x[1]*p[1]
            t3 = 4.0
            t4 = t1 + t2
            t5 = t4 + t3
            out[1] = t5
        end

        function hj(out,x,p)
            out[1] = 2*x[1]+p[1]
        end
        np = 1; nx = 1; ng = 1
        sparse_pattern = Tuple{Int64,Int64}[(1,1)]

        EAGO.build_evaluator!(low_eval, h, np, nx)
        EAGO.build_evaluator!(low_eval, h, np, nx, obj = f)
        EAGO.build_evaluator!(low_eval, h, np, nx, obj = f, constr = g, ng = 1)
        EAGO.build_evaluator!(low_eval, h, np, nx, obj = f, constr = g, ng = 1)
        EAGO.build_evaluator!(low_eval, h, np, nx, obj = f, constr = g, ng = 1, state_jac = hj)
        EAGO.build_evaluator!(low_eval, h, np, nx, obj = f, constr = g, ng = 1, state_jac = hj, user_sparse = sparse_pattern)
        @test low_eval.objective_fun == f
        @test low_eval.constraints_fun == g
        @test low_eval.state_fun == h
        @test low_eval.np == np
        @test low_eval.nx == nx
        @test low_eval.ng == ng
        @test low_eval.imp_opts.nx == nx
        @test low_eval.imp_opts.np == np
        @test low_eval.state_relax == [zero(MC{1})]
        @test low_eval.cnstr_relax == [zero(MC{1})]
        @test low_eval.state_ref_relaxation[1][1] == zero(MC{1})
        @test low_eval.state_ref_relaxation[2][1] == zero(MC{1})
        @test low_eval.state_ref_relaxation[3][1] == zero(MC{1})
        @test isnan(low_eval.last_p[1])
        @test low_eval.ref_p == [0.0]

        # initial implicit eval
        EAGO.relax_implicit!(low_eval, ytest)
        @test low_eval.obj_eval == false
        @test low_eval.cnstr_eval == false
        @test low_eval.ref_p == [7.5]
        @test isnan(low_eval.last_p[1])
        @test isapprox(low_eval.P[1].lo, 6.0, atol=1E-4) && isapprox(low_eval.P[1].hi, 9.0, atol=1E-4)
        @test isapprox(low_eval.X[1].lo, -0.78000, atol=1E-4) && isapprox(low_eval.X[1].hi, -0.4, atol=1E-4)

        test16 = low_eval.state_relax
        @test isapprox(-0.6307841683146562, low_eval.state_relax[1].cv, atol=1E-4)
        @test isapprox(-0.5209127114919797, low_eval.state_relax[1].cc, atol=1E-4)
        @test isapprox(0.09414689079323225, low_eval.state_relax[1].cv_grad[1], atol=1E-4)
        @test isapprox(0.045531201933640894, low_eval.state_relax[1].cc_grad[1], atol=1E-4)
        @test isapprox(-0.7720045045045048, low_eval.state_relax[1].Intv.lo, atol=1E-4)
        @test isapprox(-0.4, low_eval.state_relax[1].Intv.hi, atol=1E-4)

        test17 = low_eval.state_ref_relaxation
        @test isapprox(-0.728951584295428, low_eval.state_ref_relaxation[end][1].cv, atol=1E-4)
        @test isapprox(-0.4, low_eval.state_ref_relaxation[end][1].cc, atol=1E-4)
        @test isapprox(0.1470542427914973, low_eval.state_ref_relaxation[end][1].cv_grad[1], atol=1E-4)
        @test isapprox(0.0, low_eval.state_ref_relaxation[end][1].cc_grad[1], atol=1E-4)
        @test isapprox(-0.78, low_eval.state_ref_relaxation[end][1].Intv.lo, atol=1E-4)
        @test isapprox(-0.4, low_eval.state_ref_relaxation[end][1].Intv.hi, atol=1E-4)

        # second implicit eval
        EAGO.relax_implicit!(low_eval, ytest2)
        @test low_eval.obj_eval == false
        @test low_eval.cnstr_eval == false
        @test low_eval.ref_p[1] == 7.5
        @test isnan(low_eval.last_p[1])
        @test isapprox(-0.7249310591078884, low_eval.state_relax[1].cv, atol=1E-4)
        @test isapprox(-0.6653304883831199, low_eval.state_relax[1].cc, atol=1E-4)
        @test isapprox(0.09414689079323225, low_eval.state_relax[1].cv_grad[1], atol=1E-4)
        @test isapprox(0.15775515127315687, low_eval.state_relax[1].cc_grad[1], atol=1E-4)
        @test isapprox(-0.7720045045045048, low_eval.state_relax[1].Intv.lo, atol=1E-4)
        @test isapprox(-0.4, low_eval.state_relax[1].Intv.hi, atol=1E-4)
        @test isapprox(-0.728951584295428, low_eval.state_ref_relaxation[end][1].cv, atol=1E-4)
        @test isapprox(-0.4, low_eval.state_ref_relaxation[end][1].cc, atol=1E-4)
        @test isapprox(0.1470542427914973, low_eval.state_ref_relaxation[end][1].cv_grad[1], atol=1E-4)
        @test isapprox(0.0, low_eval.state_ref_relaxation[end][1].cc_grad[1], atol=1E-4)
        @test isapprox(-0.78, low_eval.state_ref_relaxation[end][1].Intv.lo, atol=1E-4)
        @test isapprox(-0.4, low_eval.state_ref_relaxation[end][1].Intv.hi, atol=1E-4)

        # test objective eval
        EAGO.relax_objective!(low_eval)
        test31 = low_eval.obj_relax
        @test low_eval.obj_eval
        @test isapprox(-0.7249310591078884, low_eval.obj_relax.cv, atol=1E-4)
        @test isapprox(-0.6653304883831199, low_eval.obj_relax.cc, atol=1E-4)
        @test isapprox(0.09414689079323225, low_eval.obj_relax.cv_grad[1], atol=1E-4)
        @test isapprox(0.15775515127315687, low_eval.obj_relax.cc_grad[1], atol=1E-4)
        @test isapprox(-0.7720045045045048, low_eval.obj_relax.Intv.lo, atol=1E-4)
        @test isapprox(-0.4, low_eval.obj_relax.Intv.hi, atol=1E-4)

        # test constraint eval
        EAGO.relax_constraints!(low_eval)
        @test isapprox(-0.7249310591078884, low_eval.cnstr_relax[1].cv, atol=1E-4)
        @test isapprox(-0.6653304883831199, low_eval.cnstr_relax[1].cc, atol=1E-4)
        @test isapprox(0.09414689079323225, low_eval.cnstr_relax[1].cv_grad[1], atol=1E-4)
        @test isapprox(0.15775515127315687, low_eval.cnstr_relax[1].cc_grad[1], atol=1E-4)
        @test isapprox(-0.7720045045045048, low_eval.cnstr_relax[1].Intv.lo, atol=1E-4)
        @test isapprox(-0.4, low_eval.cnstr_relax[1].Intv.hi, atol=1E-4)
        @test low_eval.cnstr_eval

        objval = MOI.eval_objective(low_eval, ytest)
        @test isapprox(-0.6307841683146562, objval, atol=1E-4)

        gval = zeros(1)
        MOI.eval_constraint(low_eval, gval, ytest)
        @test isapprox(-0.6307841683146562, gval[1], atol=1E-4)

        df = zeros(1)
        MOI.eval_objective_gradient(low_eval, df, ytest)
        @test isapprox(0.09414689079323225, df[1], atol=1E-4)

        jacobian_sparsity = MOI.jacobian_structure(low_eval)
        @test jacobian_sparsity[1][1] == 1
        @test jacobian_sparsity[1][2] == 1
        @test low_eval.jacobian_sparsity[1][1] == 1
        @test low_eval.jacobian_sparsity[1][2] == 1

        # throws errors
        @test_throws ErrorException MOI.hessian_lagrangian_structure(low_eval)
        @test_throws ErrorException EAGO._hessian_lagrangian_structure(low_eval)

        # no error
        @test_nowarn MOI.initialize(low_eval, [:Grad,:Jac])

        # throws error
        @test_throws ErrorException MOI.objective_expr(low_eval)
        @test_throws ErrorException MOI.constraint_expr(low_eval)

        jac_val = zeros(1,1)
        MOI.eval_constraint_jacobian(low_eval, jac_val, ytest)
        @test isapprox(0.09414689079323225, jac_val[1], atol=1E-4)

        w = 0.5
        jac_prod_val = zeros(1)
        MOI.eval_constraint_jacobian_product(low_eval, jac_prod_val, ytest, w)
        @test isapprox(jac_prod_val[1], 0.04707344539661613, atol=1E-4)

        jact_prod_val = zeros(1)
        MOI.eval_constraint_jacobian_transpose_product(low_eval, jact_prod_val, ytest, w)
        @test isapprox(jact_prod_val[1], 0.04707344539661613, atol=1E-4)

        features = MOI.features_available(low_eval)
        @test features[1] == :Grad
        @test features[2] == :Jac
    end

    # FIX ME BY ADDING IN PLACE STORAGE TO TESTS
    @testset "Upper Interval Newton" begin
        function h!(out, x,p)
            out[1] = x[1]^2 + p[1]*x[1] + 4.0
            out[2] = x[2]^2 + p[2]*x[2] + 4.0
         end
        function hj!(out,x,p)
            out[1,1] = 2.0*x[1]+ p[1]
            out[1,2] = zero(p[1])
            out[2,1] = zero(p[1])
            out[2,2] = 2.0*x[2] + p[2]
        end

        X = [EAGO.IntervalType(5.0,10.0), EAGO.IntervalType(5.0,10.0)]
        Xin = [EAGO.IntervalType(6.0,9.0), EAGO.IntervalType(6.0,9.0)]
        P = [EAGO.IntervalType(0.0,10.0), EAGO.IntervalType(0.0,10.0)]

        opts =  EAGO.interval_newton_opt()

        Y1 = zeros(Float64,2,2); J1 = zeros(IntervalType,2,2)
        Y2 = zeros(Float64,2,2); J2 = zeros(IntervalType,2,2)
        hj!(J1,X,P)
        hj!(J2,X,P)
        EAGO.preconditioner_1!(J1,Y1)
        EAGO.preconditioner_n!(J2,Y2)

        @test Y1[1] == 0.05
        @test Y2 == [0.05 0.0
                    0.0  0.05]

        bool1a = EAGO.strict_x_in_y(Xin,X,2)
        bool1b = EAGO.strict_x_in_y(X,P,2)

        bool2a = EAGO.strict_x_in_y(Xin[1],X[1])
        bool2b = EAGO.strict_x_in_y(X[1],P[1])

        bool3a = EAGO.is_equal(X,X,1E-3,2)
        bool3b = EAGO.is_equal(X,P,1E-3,2)

        @test bool1a == true
        @test bool1b == false

        @test bool2a == true
        @test bool2b == false

        @test bool3a == true
        @test bool3b == false

        function h1(out,x,p)
            out[1] = x[1]^2 + p[1]*x[1]
        end
        function hj1(out,x,p)
            out[1] = 2.0*x[1]+ p[1]
        end
        opts1 =  EAGO.interval_newton_opt(1,h1,hj1)
        X1 = [EAGO.IntervalType(5.0,10.0)]
        P1 = [EAGO.IntervalType(0.0,10.0)]

        H1 = zeros(IntervalType,1)
        J1 = zeros(IntervalType,2,2)
        Y1 = zeros(IntervalType,2,2)
        B1 = zeros(IntervalType,1)
        M1 = zeros(IntervalType,2,2)
        x_mid1 = zeros(Float64,1)
        N1 = zero(IntervalType)
        EAGO.precondition_block_1!(H1,J1,Y1,B1,M1,x_mid1,X1,P1,opts1)

        @test Y1[1] == 0.05

        function h1(out,x,p)
            out[1] = x[1]^2 + p[1]*x[1] + 4.0
        end
        function hj1(out,x,p)
            out[1] = 2.0*x[1]+ p[1]
        end
        opts1 =  EAGO.interval_newton_opt(1,h1,hj1)
        X1 = [EAGO.IntervalType(-10.0,-5.0)]
        P1 = [EAGO.IntervalType(6.0,9.0)]

        H1 = zeros(IntervalType,1)
        J1 = zeros(IntervalType,2,2)
        Y1 = zeros(Float64,2,2)
        B1 = zeros(IntervalType,1)
        M1 = zeros(IntervalType,2,2)
        x_mid1 = mid.(X1)
        N1 = zeros(IntervalType,1)
        EAGO.precondition_block_1!(H1,J1,Y1,B1,M1,x_mid1,X1,P1,opts1)
        bool4 = EAGO.interval_newton_1_kernel_nb!(M1,B1,x_mid1,X1,N1)
        EAGO.interval_newton_1_nb!(H1,J1,Y1,B1,M1,X1,P1,opts1)

        @test isapprox(H1[1].lo,  -7.25, atol=1E-4)
        @test isapprox(H1[1].hi, 15.25, atol=1E-4)
        @test isapprox(J1[1].lo, -14, atol=1E-4)
        @test isapprox(J1[1].hi, -1, atol=1E-4)
        @test isapprox(Y1[1], -0.13333333333333333, atol=1E-4)
        @test isapprox(B1[1].lo, -2.0334, atol=1E-4)
        @test isapprox(B1[1].hi, 0.966667, atol=1E-4)
        @test isapprox(M1[1].lo, 0.133333, atol=1E-4)
        @test isapprox(M1[1].hi, 1.86667, atol=1E-4)
        @test isapprox(x_mid1[1], -7.5, atol=1E-4)
        @test isapprox(N1[1].lo, -14.7501, atol=1E-4)
        @test isapprox(N1[1].hi, 7.75001, atol=1E-4)
        #=

        Y4 = EAGO.precondition_block_n(H,J,Y,B,M,x_mid,X,P,opts::interval_newton_opt)
        EAGO.interval_newton_n_kernel_nb!(M,B,x_mid,X,N,nx)
        EAGO.interval_newton_n!(x::Vector{IntervalType},p::Vector{IntervalType},opts::interval_newton_opt)
        =#

        function h(out,x,p)
            out[1] = x[1]^2 + p[1]*x[1] + 4.0
            out[2] = x[2]^2 + p[2]*x[2] + 4.0
        end
        function hj(out,x,p)
            out[1,1] = 2.0*x[1]+ p[1]
            out[1,2] = zero(p[1])
            out[2,1] = zero(p[1])
            out[2,2] = 2.0*x[2] + p[2]
        end
        X2 = [EAGO.IntervalType(-0.78,-0.4), EAGO.IntervalType(-0.78,-0.4)]
        P2 = [EAGO.IntervalType(6.0,9.0), EAGO.IntervalType(6.0,9.0)]
        opts2 = EAGO.interval_newton_opt(2,h,hj)

        H2 = zeros(IntervalType,2)
        J2 = zeros(IntervalType,2,2)
        Y2 = zeros(Float64,2,2)
        B2 = zeros(IntervalType,2)
        M2 = zeros(IntervalType,2,2)
        x_mid2 = mid.(X2)
        N2 = zeros(IntervalType,2)
        EAGO.precondition_block_n!(H2,J2,Y2,B2,M2,x_mid2,X2,P2,opts2)
        EAGO.interval_newton_n_kernel_nb!(M2,B2,x_mid2,X2,N2,2)
        EAGO.interval_newton_n!(H2,J2,Y2,B2,M2,X2,P2,opts2)

        function check_interval_array(x,lo,hi,atol)
            flag = true
            for (j,intv) in enumerate(x)
                if ~(isapprox(intv.lo,lo[j],atol=atol) && isapprox(intv.hi,hi[j],atol=atol))
                    flag = false
                    break
                end
            end
            return flag
        end

        function check_float_array(x,val,atol)
            flag = true
            for (j,intv) in enumerate(x)
                if ~isapprox(intv,val[j],atol=atol)
                    flag = false
                    break
                end
            end
            return flag
        end

        flag1 = check_interval_array(H2, [-0.929303; -0.929303], [0.828199; 0.828199], 1E-3)
        flag2 = check_interval_array(J2, [4.45666 0.0; 0.0 4.45666], [8.20001 0.0; 0.0 8.20001],1E-3)
        flag3 = check_float_array(Y2, [0.15802 0.0; 0.0 0.15802], 1E-3)

        flag4 = check_interval_array(B2, [-0.146848; -0.146848], [0.130872; 0.130872], 1E-3)
        flag5 = check_interval_array(M2, [0.70424 0.0; 0.0 0.70424], [1.29576 0.0; 0.0 1.29576], 1E-3)
        flag6 = check_float_array(x_mid2, [-0.59; -0.59], 1E-3)
        flag7 = check_interval_array(N2, [-0.772005; -0.772005], [-0.373355; -0.373355], 1E-3)
        flag8 = check_interval_array(X2, [-0.772005; -0.772005], [-0.4; -0.4], 1E-3)

        @test flag1
        @test flag2
        @test flag3
        @test flag4
        @test flag5
        @test flag6
        @test flag7
        @test flag8
    end

    @testset "Implicit Upper MidPoint Evaluator" begin

        new_node = EAGO.NodeBB(Float64[-10.0, 6.0],
                               Float64[-5.0, 9.0],
                               -Inf,Inf,0,-1,false)
        ptest = [7.5]; ptest2 = [6.5]

        function h(out,x,p)
            out[1] = x[1]^2 + p[1]*x[1] + 4.0
        end
        function hj(out,x,p)
            out[1] = 2.0*x[1]+ p[1]
        end
        X1 = [EAGO.IntervalType(-10.0,-5.0)]
        P1 = [EAGO.IntervalType(6.0,9.0)]

        f(x,p) = x[1]^2 + p[1]^2
        g(x,p) = [p[1]-7.0]

        np = 1; nx = 1; ng = 1
        y = [-7.5,7.5]
        out1 = zeros(EAGO.IntervalType,2)
        out2 = zeros(EAGO.IntervalType,1)
        mup = MidPointUpperEvaluator()
        mup1 = MidPointUpperEvaluator()

        EAGO.set_current_node!(mup,new_node)
        EAGO.set_intervals!(mup)

        EAGO.implicit_reform_5!(out1,mid.(X1),ptest,f,g,ng)
        EAGO.implicit_reform_6!(out2,mid.(X1),ptest,f)
        EAGO.build_evaluator!(mup,f,h,nx,np,hj=hj)
        EAGO.set_current_node!(mup,new_node)
        EAGO.set_intervals!(mup)
        EAGO.calc_functions!(mup,y)
        fout1 = MOI.eval_objective(mup, y)
        gout1 = Float64[]
        MOI.eval_constraint(mup, gout1, y)

        EAGO.build_evaluator!(mup1,f,h,nx,np,g=g,ng=ng,hj=hj)
        EAGO.set_current_node!(mup1,new_node)
        EAGO.set_intervals!(mup1)
        EAGO.calc_functions!(mup1,y)
        fout2 = MOI.eval_objective(mup1, y)
        gout2 = zeros(Float64,1)
        MOI.eval_constraint(mup1, gout2, y)

        flag1 = fout1 == 156.25
        flag2 = fout2 == 156.25
        flag3 = length(gout1) == 0
        flag4 = gout2[1] == 0.5
    end

    @testset "Implicit Upper NLP Evaluator" begin

        upper_eval =  EAGO.ImplicitUpperEvaluator()
        new_node = EAGO.NodeBB(Float64[5.0, 5.0, 0.0, 0.0],
                               Float64[10.0, 10.0, 10.0, 10.0],
                               -Inf,Inf,0,-1,false)
        ptest = [7.5]; ptest2 = [6.5]
        EAGO.set_current_node!(upper_eval, new_node)

        f(x,p) = x[1]^2 + x[2]^2 + p[1]^2
        g(x,p) = [x[1]
                  x[2]]
        function h(out,x,p)
            out[1] = x[1]^2 + p[1]*x[1] + 4.0
            out[2] = x[2]^2 + p[2]*x[2] + 4.0
        end
        np = 2; nx = 2; ng = 2
        y = [6.0, 7.1, 3.2, 4.3]

        out1 = zeros(7); EAGO.implicit_reform_1!(out1,y,f,g,h,np,ng,nx)
        out2 = zeros(6); EAGO.implicit_reform_2!(out2,y,g,h,np,ng,nx)
        out3 = zeros(5); EAGO.implicit_reform_3!(out3,y,f,h,np,nx)
        out4 = zeros(4); EAGO.implicit_reform_4!(out4,y,h,np,nx)

        @test isapprox(out1[1], 96.65, atol=1E-3)
        @test isapprox(out1[2], 6.0, atol=1E-3)
        @test isapprox(out1[3], 7.10, atol=1E-3)
        @test isapprox(out1[4], 59.2, atol=1E-3)
        @test isapprox(out1[5], 84.94, atol=1E-3)
        @test isapprox(out1[6], -59.2, atol=1E-3)
        @test isapprox(out1[7], -84.94, atol=1E-3)

        @test isapprox(out2[1], 6.0, atol=1E-3)
        @test isapprox(out2[2], 7.10, atol=1E-3)
        @test isapprox(out2[3], 59.2, atol=1E-3)
        @test isapprox(out2[4], 84.94, atol=1E-3)
        @test isapprox(out2[5], -59.2, atol=1E-3)
        @test isapprox(out2[6], -84.94, atol=1E-3)

        @test isapprox(out3[1], 96.65, atol=1E-3)
        @test isapprox(out3[2], 59.2, atol=1E-3)
        @test isapprox(out3[3], 84.94, atol=1E-3)
        @test isapprox(out3[4], -59.2, atol=1E-3)
        @test isapprox(out3[5], -84.94, atol=1E-3)

        @test isapprox(out4[1], 59.2, atol=1E-3)
        @test isapprox(out4[2], 84.94, atol=1E-3)
        @test isapprox(out4[3], -59.2, atol=1E-3)
        @test isapprox(out4[4], -84.94, atol=1E-3)

        EAGO.build_evaluator!(upper_eval, f, h, np, nx)
        @test upper_eval.nx == 2
        @test upper_eval.ny == 4
        @test upper_eval.ng == 0
        @test upper_eval.np == 2
        test_5e = upper_eval.value_storage
        test_6e = upper_eval.diff_result
        test_7e = upper_eval.last_y
        @test upper_eval.has_nlobj

        EAGO.calc_functions!(upper_eval, y)
        test_9e = upper_eval.diff_result
        test_10e = upper_eval.value_storage
        test_11e = upper_eval.last_y

        @test isapprox(test_9e[1,1], 12.0, atol=1E-3)
        @test isapprox(test_9e[1,2], 14.2, atol=1E-3)
        @test isapprox(test_9e[1,3], 6.4, atol=1E-3)
        @test isapprox(test_9e[2,1], 15.2, atol=1E-3)
        @test isapprox(test_9e[2,3], 6.0, atol=1E-3)
        @test isapprox(test_9e[3,2], 18.5, atol=1E-3)
        @test isapprox(test_9e[3,4], 7.1, atol=1E-3)
        @test isapprox(test_9e[4,1], -15.2, atol=1E-3)
        @test isapprox(test_9e[4,3], -6.0, atol=1E-3)
        @test isapprox(test_9e[5,2], -18.5, atol=1E-3)
        @test isapprox(test_9e[5,4], -7.1, atol=1E-3)
        @test isapprox(test_10e[1], 96.65, atol=1E-3)
        @test isapprox(test_10e[2], 59.2, atol=1E-3)
        @test isapprox(test_10e[3], 84.94, atol=1E-3)
        @test isapprox(test_10e[4], -59.2, atol=1E-3)
        @test isapprox(test_10e[5], -84.94, atol=1E-3)
        @test isapprox(test_11e[1], 6.0, atol=1E-3)
        @test isapprox(test_11e[2], 7.1, atol=1E-3)
        @test isapprox(test_11e[3], 3.2, atol=1E-3)
        @test isapprox(test_11e[4], 4.3, atol=1E-3)


        EAGO.build_evaluator!(upper_eval, f, h, np, nx; g = g, ng = 2)
        @test upper_eval.nx == 2
        @test upper_eval.ny == 4
        @test upper_eval.ng == 2
        @test upper_eval.np == 2
        test_5f = upper_eval.value_storage
        test_6f = upper_eval.diff_result
        test_7f = upper_eval.last_y
        @test upper_eval.has_nlobj

        EAGO.calc_functions!(upper_eval, y)
        test_9f = upper_eval.diff_result
        test_10f = upper_eval.value_storage

        @test isapprox(test_9f[1,1], 12.0, atol=1E-3)
        @test isapprox(test_9f[1,2], 14.2, atol=1E-3)
        @test isapprox(test_9f[1,3], 6.4, atol=1E-3)
        @test isapprox(test_9f[2,1], 1.0, atol=1E-3)
        @test isapprox(test_9f[3,2], 1.0, atol=1E-3)
        @test isapprox(test_9f[4,1], 15.2, atol=1E-3)
        @test isapprox(test_9f[4,3], 6.0, atol=1E-3)
        @test isapprox(test_9f[5,2], 18.5, atol=1E-3)
        @test isapprox(test_9f[5,4], 7.1, atol=1E-3)
        @test isapprox(test_9f[6,1], -15.2, atol=1E-3)
        @test isapprox(test_9f[6,3], -6.0, atol=1E-3)
        @test isapprox(test_9f[7,2], -18.5, atol=1E-3)
        @test isapprox(test_9f[7,4], -7.1, atol=1E-3)
        @test isapprox(test_10f[1], 96.65, atol=1E-3)
        @test isapprox(test_10f[2], 6.0, atol=1E-3)
        @test isapprox(test_10f[3], 7.1, atol=1E-3)
        @test isapprox(test_10f[4], 59.2, atol=1E-3)
        @test isapprox(test_10f[5], 84.94, atol=1E-3)
        @test isapprox(test_10f[6], -59.2, atol=1E-3)
        @test isapprox(test_10f[7], -84.94, atol=1E-3)

        fval = MOI.eval_objective(upper_eval, y)
        @test isapprox(fval, 96.65, atol=1E-3)

        gval = zeros(6)
        MOI.eval_constraint(upper_eval, gval, y)
        @test isapprox(gval[1], 6.0, atol=1E-3)
        @test isapprox(gval[2], 7.10, atol=1E-3)
        @test isapprox(gval[3], 59.2, atol=1E-3)
        @test isapprox(gval[4], 84.94, atol=1E-3)
        @test isapprox(gval[5], -59.2, atol=1E-3)
        @test isapprox(gval[6], -84.94, atol=1E-3)

        dfval = zeros(4)
        MOI.eval_objective_gradient(upper_eval, dfval, y)
        @test isapprox(dfval[1], 12.0, atol=1E-3)
        @test isapprox(dfval[2], 14.2, atol=1E-3)
        @test isapprox(dfval[3], 6.4, atol=1E-3)
        @test isapprox(dfval[4], 0.0, atol=1E-3)

        jac_struct = MOI.jacobian_structure(upper_eval)
        @test jac_struct[1][1] == 1
        @test jac_struct[1][2] == 1
        @test jac_struct[2][1] == 1
        @test jac_struct[2][2] == 2
        @test jac_struct[3][1] == 1
        @test jac_struct[3][2] == 3
        @test jac_struct[4][1] == 1
        @test jac_struct[4][2] == 4

        feat_sym = MOI.features_available(upper_eval)
        @test feat_sym[1] == :Grad
        @test feat_sym[2] == :Jac

        dg = zeros(24)
        MOI.eval_constraint_jacobian(upper_eval, dg, y)
        #=
        @test dg ==  [1.0  0.0   0.0   0.0
                      0.0  1.0   0.0   0.0
                      15.2  0.0  6.0  0.0
                      0.0  18.5  0.0  7.1
                     -15.2  0.0 -6.0  0.0
                      0.0 -18.5  0.0 -7.1]
        =#
        #out = zeros(6); w = fill(0.5,(6,))
        #MOI.eval_constraint_jacobian_product(upper_eval, out, y, w)
        #MOI.eval_constraint_jacobian_transpose_product(upper_eval, y, p, w)

        # should error
        @test_throws ErrorException MOI.hessian_lagrangian_structure(upper_eval)
        @test_throws ErrorException MOI.objective_expr(upper_eval)
        @test_throws ErrorException MOI.constraint_expr(upper_eval)
        #@test_throws ErrorException _hessian_lagrangian_structure(upper_eval)
    end
    @testset "Implicit Optimizer Full Routine" begin
        using MathOptInterface, EAGO

        f(x,p) = x[1]

        function h(out,x,p)
            out[1] = x[1]^2 + x[1]*p[1] + 4.0
        end
        function hj(out,x,p)
            out[1] = 2.0*x[1]+p[1]
        end

        pl = [6.0]; pu = [9.0]; xl = [-0.78]; xu = [-0.4]
        opt = EAGO.Optimizer(absolute_tolerance=1e-3)
        var, opt = solve_implicit(f, h, xl, xu, pl, pu, opt, hj, nothing)
        var_value_2 = MOI.get(opt, MOI.VariablePrimal(), var[2])
        tstatus = MOI.get(opt, MOI.TerminationStatus())
        pstatus = MOI.get(opt, MOI.PrimalStatus())
        objval = MOI.get(opt, MOI.ObjectiveValue())

        @test isapprox(var_value_2, 6.00, atol=1e-1)
        @test isapprox(objval, -0.78, atol=1e-1)
    end
end
#=
#using MathOptInterface
#const MOI = MathOptInterface
using MathOptInterface, EAGO
const MOI = MathOptInterface
f(x,p) = x[1]
h(x,p) = [x[1]^2 + x[1]*p[1] + 4.0]
hj(x,p) = [2.0*x[1]+p[1]]

pl = [6.0]; pu = [9.0]; xl = [-10.0]; xu = [-5.0]
opt = EAGO.Optimizer(AbsoluteTolerance=1e-1)
var, opt = solve_implicit(f, h, xl, xu, pl, pu, opt, hj, nothing, upper = :MidPoint)
var_value_1 = MOI.get(opt, MOI.VariablePrimal(), var[1])
var_value_2 = MOI.get(opt, MOI.VariablePrimal(), var[2])
tstatus = MOI.get(opt, MOI.TerminationStatus())
pstatus = MOI.get(opt, MOI.PrimalStatus())
objval = MOI.get(opt, MOI.ObjectiveValue())

println("objval: $objval")
flag1 = isapprox(var_value_1, -0.78, atol=1e-1)
flag2 = isapprox(objval, -0.78, atol=1e-1)
=#
