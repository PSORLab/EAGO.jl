
#=
Basic idea behind syntax
mSIP = SIPModel()

@decision_variable(mSIP, 0.0 <= x <= 1.0)
@decision_variable(mSIP, z, Bin)
@uncertain_variable(mSIP, p, 0.0 <= x <= 1.0)


@constraint(mSIP, ...)
@objective(mSIP, Min, ...)

@NLconstraint(mSIP, ...)
@NLobjective(mSIP, Min, ...)

optimize!()
=#

using JuMP
using MathOptInterface

m = Model()
@variable(m, 0 <= y <= 1)
@variable(m, x)
@variable(m, a, Bin)
@variable(m, q[i=1:2])
@constraint(m, y^2 + y + x<= 0)
@constraint(m, [x, y-1, y-2] in SecondOrderCone())
@NLconstraint(m, sin(x) + cos(y) <= 0.0)
@constraint(m, 2x - 1 ⟂ x)
@constraint(m, q in SOS2([3,5]))
@constraint(m, a => {x + y <= 1})
@SDconstraint(m, [x 2x; 3x 4x] >= ones(2, 2))

A = [1 2; 3 4]
b = [5,6]
@constraint(m, con, A * x .== b)

list = list_of_constraint_types(m)

cons0 = all_constraints(m, VariableRef, MOI.LessThan{Float64})
cons1 = all_constraints(m, VariableRef, MOI.GreaterThan{Float64})
cons2 = all_constraints(m, GenericQuadExpr{Float64,VariableRef}, MOI.LessThan{Float64})
cons2 = all_constraints(m, GenericQuadExpr{Float64,VariableRef}, MOI.LessThan{Float64})

cons0_1 = cons0[1]
cons1_1 = cons1[1]
cons2_1 = cons2[1]
cons3_1 = cons2[1]

out = constraint_object(cons2_1).func
l_terms = linear_terms(out)
q_terms = quad_terms(out)
#all_consts = all_constraints(m, list[1][1], list[1][1])