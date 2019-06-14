# circle global lib

m = Model(with_optimizer(EAGO.Optimizer))

# ----- Variables ----- #
@variable(m, objvar)
x_Idx = Any[1, 2]
@variable(m, x[x_Idx])
JuMP.set_lower_bound(objvar, 0.0)
JuMP.set_upper_bound(objvar, 300.0)
JuMP.set_lower_bound(x[1], 0.0)
JuMP.set_lower_bound(x[2], 0.0)
JuMP.set_upper_bound(x[1], 10.0)
JuMP.set_upper_bound(x[2], 10.0)

# ----- Constraints ----- #
@constraint(m, e1,  (2.545724188-x[1])^2 + (9.983058643-x[2])^2 - (objvar)^2 <= 0.0)
@constraint(m, e2,  (8.589400372-x[1])^2 + (6.208600402-x[2])^2 - (objvar)^2 <= 0.0)
@constraint(m, e3,  (5.953378204-x[1])^2 + (9.920197351-x[2])^2 - (objvar)^2 <= 0.0)
@constraint(m, e4,  (3.710241136-x[1])^2 + (7.860254203-x[2])^2 - (objvar)^2 <= 0.0)
@constraint(m, e5,  (3.629909053-x[1])^2 + (2.176232347-x[2])^2 - (objvar)^2 <= 0.0)
@constraint(m, e6,  (3.016475803-x[1])^2 + (6.757468831-x[2])^2 - (objvar)^2 <= 0.0)
@constraint(m, e7,  (4.148474536-x[1])^2 + (2.435660776-x[2])^2 - (objvar)^2 <= 0.0)
@constraint(m, e8,  (8.706433123-x[1])^2 + (3.250724797-x[2])^2 - (objvar)^2 <= 0.0)
@constraint(m, e9,  (1.604023507-x[1])^2 + (7.020357481-x[2])^2 - (objvar)^2 <= 0.0)
@constraint(m, e10,  (5.501896021-x[1])^2 + (4.918207429-x[2])^2 - (objvar)^2 <= 0.0)

# ----- Objective ----- #
@objective(m, Min, objvar)

JuMP.optimize!(m)

fval = JuMP.objective_value(m)
status_term = JuMP.termination_status(m)
status_prim = JuMP.primal_status(m)

#=
@testset "QP Problem #2: circle (global library)" begin

    m = Model(with_optimizer(EAGO.Optimizer))

    # ----- Variables ----- #
    @variable(m, objvar)
    x_Idx = Any[1, 2]
    @variable(m, x[x_Idx])
    setlowerbound(objvar, 0.0)


    # ----- Constraints ----- #
    @NLconstraint(m, e1,  (2.545724188-x[1])^2+ (9.983058643-x[2])^2- (objvar)^2 <= 0.0)
    @NLconstraint(m, e2,  (8.589400372-x[1])^2+ (6.208600402-x[2])^2- (objvar)^2 <= 0.0)
    @NLconstraint(m, e3,  (5.953378204-x[1])^2+ (9.920197351-x[2])^2- (objvar)^2 <= 0.0)
    @NLconstraint(m, e4,  (3.710241136-x[1])^2+ (7.860254203-x[2])^2- (objvar)^2 <= 0.0)
    @NLconstraint(m, e5,  (3.629909053-x[1])^2+ (2.176232347-x[2])^2- (objvar)^2 <= 0.0)
    @NLconstraint(m, e6,  (3.016475803-x[1])^2+ (6.757468831-x[2])^2- (objvar)^2 <= 0.0)
    @NLconstraint(m, e7,  (4.148474536-x[1])^2+ (2.435660776-x[2])^2- (objvar)^2 <= 0.0)
    @NLconstraint(m, e8,  (8.706433123-x[1])^2+ (3.250724797-x[2])^2- (objvar)^2 <= 0.0)
    @NLconstraint(m, e9,  (1.604023507-x[1])^2+ (7.020357481-x[2])^2- (objvar)^2 <= 0.0)
    @NLconstraint(m, e10,  (5.501896021-x[1])^2+ (4.918207429-x[2])^2- (objvar)^2 <= 0.0)


    # ----- Objective ----- #
    @objective(m, Min, objvar)

    JuMP.optimize!(m)

    fval = JuMP.objective_value(m)
    status_term = JuMP.termination_status(m)
    status_prim = JuMP.primal_status(m)

    @test isapprox(fval,4.5742,atol=1E-3)
    @test status_term == MOI.Success
    @test status_prim == MOI.FeasiblePoint
end
=#
