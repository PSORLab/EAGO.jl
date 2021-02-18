
#=
using JuMP

function pseudo_copy(m::JuMP.Model)

    nlp_new = deepcopy(m.nlp_data)
    nlp_new.evaluator = NLPEvaluator(m)
    MOI.initialize(nlp_new.evaluator, Symbol[:ExprGraph])

    obj_expr = MOI.objective_expr(nlp_new.evaluator)
    con_expr = MOI.constraint_expr(nlp_new.evaluator, 1)

    @show obj_expr
    @show con_expr

    return nothing
end

m = Model()
@variable(m, x[1:3])
@variable(m, y)
@NLconstraint(m, x[1] + y^2 <= 0.0)
@NLobjective(m, Min, x[2]*x[3]*y^2)

nlp_data = m.nlp_data

out = pseudo_copy(m)
=#
