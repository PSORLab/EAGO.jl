module TestDAGcntr

using Compat
using Compat.Test
using IntervalArithmetic
using EAGO

single_expr = :(x[1]+x[2]^2-exp(x[1]))
tape_out = Generate_Tape(single_expr,2,-2,4,Interval{Float64})
dag_out = getDAG(single_expr,[1,2],Interval{Float64})
@test dag_out[1] == Int64[1; 2; 3; 4; 5; 6; 7]
@test dag_out[3][4] == Int64[4; 5]

#=
generates the tape list for the constraints -1 <= (x[1]+x[2]^2-exp(x[1]) <= 2
and 2 <= :(x[1]*x[2]-x[5]) <= 4 and runs the contractor 6 times in place on the
variable box X = [1,20]^5
=#
mult_expr = [:(x[1]+x[2]^2)
             :(x[1]+x[5])
            ]
X = [Interval(-10.0,20.0) for i=1:5]
ftapelist_out = Generate_Fixed_TapeList(mult_expr,4,[-8.0,0.5],[2.0,4.5],[[5.0]],Interval{Float64})
DAGContractor!(X,ftapelist_out,6)
@test -10-1E-4 <= X[1].lo <= -10+1E-4
@test 0-1E-4 <= X[2].lo <= 0+1E-4
@test -10-1E-4 <= X[3].lo <= -10+1E-4
@test -10-1E-4 <= X[4].lo <= -10+1E-4
@test 5-1E-4 <= X[5].lo <= 5+1E-4
@test 20-1E-4 <= X[1].hi <= 20+1E-4
@test 20-1E-4 <= X[2].hi <= 20+1E-4
@test 20-1E-4 <= X[3].hi <= 20+1E-4
@test 20-1E-4 <= X[4].hi <= 20+1E-4
@test 5-1E-4 <= X[5].hi <= 5+1E-4

end
