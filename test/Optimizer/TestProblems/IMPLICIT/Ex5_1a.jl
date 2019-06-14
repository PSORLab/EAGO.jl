opt1 = EAGO.Optimizer(absolute_tolerance = 0.001, relative_tolerance = 0.001)
obj1(x,p) = x[1]
function h1!(out,x,p)
    out[1] = x[1]^2 + x[1]*p[1] + 4.0
end
function hjac1!(out,x,p)
    out[1] = 2.0*x[1] + p[1]
end
g1(x,p) = [-x[1] - 6.0]
pl1 = [6.0]; pu1 = [9.0];
xl1 = [-10.0]; xu1 = [-5.0]
var1, opt1 = solve_implicit(obj1, h1!, xl1, xu1, pl1, pu1, opt1, hjac1!, g1)

pval1 = MOI.get(opt1, MOI.VariablePrimal(), var1[2])
fval1 = MOI.get(opt1, MOI.ObjectiveValue())
tstatus1 = MOI.get(opt1, MOI.TerminationStatus())
pstatus1 = MOI.get(opt1, MOI.PrimalStatus())

println("pval1: $pval1")
println("fval1: $fval1")

p_bool1 = isapprox(pval1, 9.0, atol = 0.01)
f_bool1 = isapprox(fval1, -8.53113, atol = 0.001)
t_bool1 = (tstatus1 == MOI.OPTIMAL)
ps_bool1 = (pstatus1 == MOI.FEASIBLE_POINT)

@test p_bool1
@test f_bool1
@test t_bool1
@test ps_bool1
