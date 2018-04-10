#workspace()

using IntervalArithmetic
using Clp
using EAGODomainReduction

type t_solver
    LP_solver
end

type test_opt
    numVar
    numConstr
    f
    g
    solver::t_solver
    gL
    gU
    gL_loc
    gU_loc
end

f(x) = x[1]^2+x[2]^2
g(x) = [x[1]+x[2];
        x[1]-x[2]]

in_opt = t_solver(ClpSolver())
opt = test_opt(2,2,f,g,in_opt,[2 -1],[3 0],[1,2],[1,2])
X = [Interval(0.0,2.0) for i=1:2]
STD_Linear_RR!(X,[opt],7.0)
println("Standard Range Reduced X: ", X)
