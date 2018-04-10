module Options_Test

using Compat
using Compat.Test
using EAGOBranchBound

x = BnBSolver()

set_Branch_Scheme!(x,"best")
@test x.Node_Select == EAGOBranchBound.NS_best
@test x.Branch_Sto == EAGOBranchBound.BM_depth_best!

set_Branch_Scheme!(x,"breadth")
@test x.Node_Select == EAGOBranchBound.NS_depth_breadth
@test x.Branch_Sto == EAGOBranchBound.BM_breadth!

set_Branch_Scheme!(x,"depth")
@test x.Node_Select == EAGOBranchBound.NS_depth_breadth
@test x.Branch_Sto == EAGOBranchBound.BM_depth_best!

set_Bisect_Func!(x,"relative midpoint",-1)
@test x.Bisect_Func == EAGOBranchBound.Bisect_Rel
set_Bisect_Func!(x,"absolute midpoint",-1)
@test x.Bisect_Func == EAGOBranchBound.Bisect_Abs

set_Verbosity!(x,"Normal")
@test x.Verbosity == "Normal"

@test_throws ErrorException set_Branch_Scheme!(x,"random")
@test_throws ErrorException  set_Verbosity!(x,"random")
@test_throws ErrorException  set_Bisect_Func!(x,"random",-1)
@test_throws ErrorException  set_Bisect_Func!(x,"random",1)

@test_nowarn set_Bisect_Func!(x,"relative midpoint",3)
@test_nowarn set_Bisect_Func!(x,"absolute midpoint",3)

end
