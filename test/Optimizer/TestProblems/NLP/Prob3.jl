using EAGO, JuMP
m = Model(with_optimizer(EAGO.Optimizer, cut_max_iterations = 10, verbosity = 0))

# ----- Variables ----- #
@variable(m, objvar)
x_Idx = Any[1, 2, 3]
@variable(m, x[x_Idx])
JuMP.set_lower_bound(objvar, -100.0)
JuMP.set_upper_bound(objvar, 100.0)
JuMP.set_lower_bound(x[1], -5.0)
JuMP.set_upper_bound(x[1], 5.0)
JuMP.set_lower_bound(x[2], -5.0)
JuMP.set_upper_bound(x[2], 5.0)
JuMP.set_lower_bound(x[3], -100.0)
JuMP.set_upper_bound(x[3], 100.0)


# ----- Constraints ----- #
@constraint(m, e1, -x[3]+objvar == 0.0)
@NLconstraint(m, e2, 2*(x[2])^2+4*x[1]*x[2]-42*x[1]+4*(x[1])^3-x[3] <= 14.0)
@NLconstraint(m, e3, (-2*(x[2])^2)-4*x[1]*x[2]+42*x[1]-4*(x[1])^3-x[3] <= -14.0)
@NLconstraint(m, e4, 2*(x[1])^2+4*x[1]*x[2]-26*x[2]+4*(x[2])^3-x[3] <= 22.0)
@NLconstraint(m, e5, (-2*(x[1])^2)-4*x[1]*x[2]+26*x[2]-4*(x[2])^3-x[3] <= -22.0)

# ----- Objective ----- #
@objective(m, Min, objvar)

JuMP.optimize!(m)

fval = JuMP.objective_value(m)
status_term = JuMP.termination_status(m)
status_prim = JuMP.primal_status(m)

@testset "NLP Problem #4: ex14_1_1 (global library)" begin

    m = Model(with_optimizer(EAGO.Optimizer))

    # ----- Variables ----- #
    @variable(m, objvar)
    x_Idx = Any[1, 2, 3]
    @variable(m, x[x_Idx])
    setlowerbound(x[1], -5.0)
    setupperbound(x[1], 5.0)
    setlowerbound(x[2], -5.0)
    setupperbound(x[2], 5.0)


    # ----- Constraints ----- #
    @constraint(m, e1, -x[3]+objvar == 0.0)
    @NLconstraint(m, e2, 2* (x[2])^2+4*x[1]*x[2]-42*x[1]+4* (x[1])^3-x[3] <= 14.0)
    @NLconstraint(m, e3, (-2* (x[2])^2)-4*x[1]*x[2]+42*x[1]-4* (x[1])^3-x[3] <= -14.0)
    @NLconstraint(m, e4, 2* (x[1])^2+4*x[1]*x[2]-26*x[2]+4* (x[2])^3-x[3] <= 22.0)
    @NLconstraint(m, e5, (-2* (x[1])^2)-4*x[1]*x[2]+26*x[2]-4* (x[2])^3-x[3] <= -22.0)


    # ----- Objective ----- #
    @objective(m, Min, objvar)

    fval = JuMP.objective_value(m)
    status_term = JuMP.termination_status(m)
    status_prim = JuMP.primal_status(m)

    @test isapprox(fval,0.000,atol=1E-3)
    @test status_term == MOI.Success
    @test status_prim == MOI.FeasiblePoint
end
