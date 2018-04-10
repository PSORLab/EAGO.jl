module TestRR

using Compat
using Compat.Test
using IntervalArithmetic
using Clp
using EAGO

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

@test 0.5-1E-4 <= X[1].lo <= 0.5+1E-4
@test 1.25-1E-4 <= X[1].hi <= 1.25+1E-4
@test 1.0-1E-4 <= X[2].lo <= 1.0+1E-4
@test 1.75-1E-4 <= X[2].hi <= 1.75+1E-4

end
