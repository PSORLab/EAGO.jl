m = Model(with_optimizer(EAGO.Optimizer))

# ----- Variables ----- #
x_Idx = Any[1, 2]
@variable(m, x[x_Idx])
JuMP.set_lower_bound(x[1], 0.0)
JuMP.set_lower_bound(x[2], 0.0)
JuMP.set_upper_bound(x[1], 3.0)
JuMP.set_upper_bound(x[2], 4.0)

# ----- Constraints ----- #
@NLconstraint(m, e2, x[2] - 8*(x[1])^2 + 8*(x[1])^3-2*(x[1])^4 <= 2.0)
@NLconstraint(m, e3, 32*(x[1])^3-4*(x[1])^4-88*(x[1])^2+96*x[1]+x[2] <= 36.0)


# ----- Objective ----- #
@objective(m, Min, -x[1]-x[2])

JuMP.optimize!(m)

#jackend_eval_obj = jbackend_eval.objective
#jackend_eval_obj = jbackend_eval.constraints

#fval = JuMP.objective_value(m)
#status_term = JuMP.termination_status(m)
#status_prim = JuMP.primal_status(m)

#=
@testset "NLP Problem #5: ex4_1_9 (global library)" begin

    m = Model(with_optimizer(EAGO.Optimizer))

    # ----- Variables ----- #
    @variable(m, objvar)
    x_Idx = Any[1, 2]
    @variable(m, x[x_Idx])
    setlowerbound(x[1], 0.0)
    setlowerbound(x[2], 0.0)
    setupperbound(x[1], 3.0)
    setupperbound(x[2], 4.0)
    setupperbound(objvar, -20.0)
    setupperbound(objvar, 20.0)


    # ----- Constraints ----- #
    @constraint(m, e1, x[1]+x[2]+objvar == 0.0)
    @NLconstraint(m, e2, 8* (x[1])^3-2* (x[1])^4-8* (x[1])^2+x[2] <= 2.0)
    @NLconstraint(m, e3, 32* (x[1])^3-4* (x[1])^4-88* (x[1])^2+96*x[1]+x[2] <= 36.0)


    # ----- Objective ----- #
    @objective(m, Min, objvar)

    fval = JuMP.objective_value(m)
    status_term = JuMP.termination_status(m)
    status_prim = JuMP.primal_status(m)

    @test isapprox(fval,-5.5080,atol=1E-3)
    @test status_term == MOI.Success
    @test status_prim == MOI.FeasiblePoint
end
=#
