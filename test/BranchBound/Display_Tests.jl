module BBDisplayTest

using Compat
using Compat.Test
using IntervalArithmetic
using EAGO

S = BnBSolver()
S.Verbosity = "Full"

@test_nowarn EAGO.print_node!(S,Int64(1),3.0,[Interval(1.0,2.0)])
@test_nowarn EAGO.print_results!(S,2.0,[1.5],false,true)
@test_nowarn EAGO.print_results!(S,2.0,[1.5],false,false)

end
