#workspace()

using IntervalArithmetic
using EAGODomainReduction

single_expr = :(x[1]+x[2]^2-exp(x[1]))
tape_out = Generate_Tape(single_expr,2,-2,4)
dag_out = getDAG(single_expr,[1,2])

#=
generates the tape list for the constraints -1 <= (x[1]+x[2]^2-exp(x[1]) <= 2
and 2 <= :(x[1]*x[2]-x[5]) <= 4 and runs the contractor 6 times in place on the
variable box X = [1,20]^5
=#
mult_expr = [:(x[1]+x[2]^2)
             :(x[1]+x[5])
            ]
tapelist_out = Generate_TapeList(mult_expr,5,[-1.0,2.0],[2.0,3.0])
X = [Interval(0.5,20.0) for i=1:5]
DAGContractor!(X,tapelist_out,6)
println("contracted X value: ", X)

#=
generates the tape list for the constraints -1 <= (x[1]+x[2]^2-exp(x[1]) <= 2
and 2 <= :(x[1]*x[2]-x[5]) <= 4 with x[5]==5 and runs the contractor 6 times in place on the
variable box X = [1,20]^5
=#

X = [Interval(-10.0,20.0) for i=1:5]
ftapelist_out = Generate_Fixed_TapeList(mult_expr,4,[-8.0,0.5],[2.0,4.5],[[5.0]])
DAGContractor!(X,ftapelist_out,6)
println("contracted X value: ", X)
