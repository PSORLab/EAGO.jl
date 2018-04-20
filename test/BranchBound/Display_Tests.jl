module DisplayTestBB

using Compat
using Compat.Test
using EAGO
using IntervalArithmetic

S = BnBSolver()
S.Verbosity = "Full"

println("test #1:")
@test_nowarn EAGO.print_node!(S,1,3.0,[Interval(1.0,2.0)])
println("test #2:")
@test_nowarn EAGO.print_results!(S,2.0,[1.5],false,true)
println("test #3:")
@test_nowarn EAGO.print_results!(S,2.0,[1.5],false,false)

end
