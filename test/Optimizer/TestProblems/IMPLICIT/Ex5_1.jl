#=
@testset "Begin Implicit Solver 5.1 Example" begin
end
=#

opt = EAGO.Optimizer(absolute_tolerance = 0.001, relative_tolerance = 0.001)
obj(x,p) = x[1]
function h!(out,x,p)
    out[1] = x[1]^2 + x[1]*p[1] + 4.0
end
function hjac!(out,x,p)
    out[1] = 2.0*x[1] + p[1]
end
pl = [6.0]; pu = [9.0];
xl = [-10.0]; xu = [-5.0]
var, opt = solve_implicit(obj, h!, xl, xu, pl, pu, opt, hjac!, nothing)

pval = MOI.get(opt, MOI.VariablePrimal(), var[2])
fval = MOI.get(opt, MOI.ObjectiveValue())
tstatus = MOI.get(opt, MOI.TerminationStatus())
pstatus = MOI.get(opt, MOI.PrimalStatus())

println("pval: $pval")
println("fval: $fval")

p_bool = isapprox(pval, 8.808071047067642, atol = 0.01)
f_bool = isapprox(fval, -8.327749184410186, atol = 0.001)
t_bool = (tstatus == MOI.OPTIMAL)
ps_bool = (pstatus == MOI.FEASIBLE_POINT)

@test p_bool
@test f_bool
@test t_bool
@test ps_bool
