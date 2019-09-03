@testset "NLP Evaluator" begin

        m = Model(with_optimizer(EAGO.Optimizer))
        @variable(m, x)
        @variable(m, y)

        @NLobjective(m, Min, exp(1 - x) ^ 2 + 100 * (y - x ^ 2) ^ 2)

        @constraint(m, x^2 + y <= 10)
        @constraint(m, x + y == 10)
        @constraint(m, y >= 0)

        @NLconstraint(m, log(y - x ^ 2) <= 0)

        source_evaluator = JuMP.NLPEvaluator(m)
        MOI.initialize(source_evaluator , Symbol[:Grad])

        opt = m.moi_backend.optimizer.model.optimizer
        built_evaluator = EAGO.build_nlp_evaluator(2, source_evaluator, opt, true)

        # Add current node and define point
        built_evaluator.current_node = EAGO.NodeBB(Float64[1.0,5.0], Float64[2.0,6.0], -Inf, Inf, 2, 1, true)
        xpoint = Float64[1.5,5.5]

        built_evaluator.constraints_lbd = [-Inf]
        built_evaluator.constraints_ubd = [0.0];

        user_operators = built_evaluator.m.nlp_data.user_operators
        nlp_data = built_evaluator.m.nlp_data
        user_input_buffer = built_evaluator.jac_storage

        subgrad_tighten = false
        first_eval_flag = false

        EAGO.forward_eval(built_evaluator.objective.setstorage, built_evaluator.objective.numberstorage,
                          built_evaluator.objective.numvalued,
                          built_evaluator.objective.nd, built_evaluator.objective.adj,
                          built_evaluator.objective.const_values, built_evaluator.parameter_values,
                          built_evaluator.current_node, xpoint, built_evaluator.subexpression_values_flt,
                          built_evaluator.subexpression_values_set, built_evaluator.subexpression_isnum,
                          user_input_buffer, subgrad_tighten, built_evaluator.objective.tpdict,
                          built_evaluator.objective.tp1storage, built_evaluator.objective.tp2storage,
                          built_evaluator.objective.tp3storage, built_evaluator.objective.tp4storage,
                          first_eval_flag, user_operators = user_operators)

        EAGO.forward_eval_all(built_evaluator,xpoint)

        EAGO.reverse_eval(built_evaluator.objective.setstorage,
                          built_evaluator.objective.numberstorage,
                          built_evaluator.objective.numvalued,
                          built_evaluator.subexpression_isnum,
                          built_evaluator.subexpression_values_set,
                          built_evaluator.objective.nd,
                          built_evaluator.objective.adj, xpoint,
                          built_evaluator.current_node, subgrad_tighten)

        EAGO.reverse_eval_all(built_evaluator,xpoint)

        EAGO.forward_reverse_pass(built_evaluator,xpoint)

        xstar1 = xpoint
        gstar1 = [1.2]
        dfstar1 = xpoint
        Jstar1 = dfstar1'

        temp1 = MOI.features_available(built_evaluator)
        temp2 = MOI.eval_objective(built_evaluator, xstar1)
        temp3 = MOI.eval_constraint(built_evaluator, gstar1, xstar1)
        temp4 = MOI.eval_objective_gradient(built_evaluator, dfstar1, xstar1)
        temp5 = MOI.jacobian_structure(built_evaluator)
        temp6 = MOI.eval_constraint_jacobian(built_evaluator, Jstar1, xstar1)
        temp7 = MOI.eval_constraint_jacobian_product(built_evaluator, y, xstar1, xpoint)
        temp8 = MOI.eval_constraint_jacobian_transpose_product(built_evaluator, y, xstar1, xpoint)

        nlconstraint = built_evaluator

        @test temp1[1] == :Grad
        @test temp1[2] == :Jac
        @test isapprox(dfstar1[1],-5099.264424,atol=1E-3)  # FAIL
        @test dfstar1[2] == 3400.0                         # FAIL
        @test temp5[1][1] == 1
        @test temp5[1][2] == 1
        @test temp5[2][1] == 1
        @test temp5[2][2] == 2
end
