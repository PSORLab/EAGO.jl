m = Model(with_optimizer(EAGO.Optimizer, InitialRelaxedOptimizer = CPLEX.Optimizer()))

# ----- Variables ----- #
@variable(m, objvar)
x_Idx = Any[1, 2, 3, 4, 5]
@variable(m, x[x_Idx])
set_lower_bound(x[5], 0.0)
set_lower_bound(x[1], 0.496094)
#set_lower_bound(x[1], 0.0)
set_lower_bound(x[4], 0.0)
set_lower_bound(x[2], 0.0)
set_lower_bound(x[3], 0.0)
set_upper_bound(x[1], 1.0)
set_upper_bound(x[2], 1.0)
set_upper_bound(x[3], 1.0)
set_upper_bound(x[4], 1.0)
set_upper_bound(x[5], 1.0)
set_lower_bound(objvar, -30.0)
set_upper_bound(objvar, -0.234375)
#set_upper_bound(objvar, 30.0)

# ----- Constraints ----- #
@constraint(m, e1, -(42*x[1]+44*x[2]+45*x[3]+47*x[4]+47.5*x[5]-0.5*(100*x[1]*x[1]+100*x[2]*x[2]+100*x[3]*x[3]+100*x[4]*x[4]+100*x[5]*x[5]))+objvar == 0.0)
@constraint(m, e2, 20*x[1]+12*x[2]+11*x[3]+7*x[4]+4*x[5] <= 40.0)

# ----- Objective ----- #
@objective(m, Min, objvar)

JuMP.optimize!(m)

fval = JuMP.objective_value(m)
status_term = JuMP.termination_status(m)
status_prim = JuMP.primal_status(m)

#=
@testset "QP Problem #3: ex2_1_1 (global library)" begin

    m = Model(with_optimizer(EAGO.Optimizer))

    # ----- Variables ----- #
    @variable(m, objvar)
    x_Idx = Any[1, 2, 3, 4, 5]
    @variable(m, x[x_Idx])
    setlowerbound(x[5], 0.0)
    setlowerbound(x[1], 0.0)
    setlowerbound(x[4], 0.0)
    setlowerbound(x[2], 0.0)
    setlowerbound(x[3], 0.0)
    setupperbound(x[1], 1.0)
    setupperbound(x[2], 1.0)
    setupperbound(x[3], 1.0)
    setupperbound(x[4], 1.0)
    setupperbound(x[5], 1.0)

    # ----- Constraints ----- #
    @NLconstraint(m, e1, -(42*x[1]-0.5*(100*x[1]*x[1]+100*x[2]*x[2]+100*x[3]*x[3]+100*x[4]*x[4]+100*x[5]*x[5])+44*x[2]+45*x[3]+47*x[4]+47.5*x[5])+objvar == 0.0)
    @constraint(m, e2, 20*x[1]+12*x[2]+11*x[3]+7*x[4]+4*x[5] <= 40.0)

    # ----- Objective ----- #
    @objective(m, Min, objvar)

    fval = JuMP.objective_value(m)
    status_term = JuMP.termination_status(m)
    status_prim = JuMP.primal_status(m)

    @test isapprox(fval,-17.0,atol=1E-3)
    @test status_term == MOI.Success
    @test status_prim == MOI.FeasiblePoint
end
=#
