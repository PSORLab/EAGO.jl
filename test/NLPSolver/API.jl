module Test_Explicit_NLP

using Compat
using Compat.Test
using EAGO

@testset " Explicit API Functionality" begin
# Usage or any other factoriable function (currently doesn't support x in if/while)
# Define pretty much any function that the library supports (error will flag if unsupported)
ffunc(x) = x[1]^3+x[2]           # objective function f(x)
ggfunc(x) = [x[1], -x[2]]        # constraint function g(x,p) <= 0
hhfunc(x) = [x[1]-2*x[2]]        # constraint function h(x,p) == 0

# Sets up and solves unconstrained explicit problem (Works)
f1, x1, feas1 = Optimize_Script(ffunc, Float64[1.0,2.0], Float64[2.0,3.0])
@test feas1 == :Optimal
@test isapprox(x1[1],1.0,atol=1e-5)
@test isapprox(x1[2],2.0,atol=1e-5)
@test isapprox(f1,3.0,atol=1e-5)

# Sets up and solves constrained explicit problem
f2, x2, feas2 = Optimize_Script(ffunc, Float64[1.0,2.0], Float64[2.0,3.0],
                                g = ggfunc, h = hhfunc)
@test feas2 == :Infeasible

pffunc(x,p) = x[1]^3+x[2]+p[1]*p[2]*x[1]*x[2]      # objective function f(x,p)
pgfunc(x,p) = [1.5-x[1], p[2]*x[1]-x[2]]              # inequality constraint function g(x,p) <= 0

# Sets up and solves constained explicit problem with parameter values (ERROR SHOULD BE FEASIBLE)
f3, x3, feas3 = Optimize_Script(pffunc, Float64[1.0,2.0], Float64[2.0,3.0],
                                g = pgfunc, p = [2.0, 1.0])

@test feas3 == :Optimal
@test isapprox(x3[1],1.5,atol=1e-5)
@test isapprox(x3[2],2.0,atol=1e-5)
@test isapprox(f3,11.374999,atol=1e-5)

end
end
