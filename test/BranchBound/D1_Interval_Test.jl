module D1_Interval_Test

using Compat
using Compat.Test
using IntervalArithmetic
using EAGO

X = [Interval(1,5)]
model = BnBModel(X)
solver = BnBSolver()
set_to_default!(solver)
solver.Verbosity = "Normal"
solver.opt = [x -> x[1]^2]
solver.Lower_Prob = (X,k,p,opt,temp) -> [opt[1](X).lo,mid.(X),true,[]]
solver.Upper_Prob = (X,k,p,opt,temp) -> [opt[1](X).hi,mid.(X),true,[]]
a = solver.Lower_Prob(X,1,1,solver.opt,[])
b = solver.Upper_Prob(X,1,1,solver.opt,[])
solveBnB!(solver,model)

@test -1E-4 <= getobjval(model)-1 <= 1E-4

solver.iter_lim = true
solver.max_iter = 10
solver.Lower_Prob = (X,k,p,opt,temp) -> [opt[1](X).lo,mid.(X),false,[]]
model2 = BnBModel(X)
solveBnB!(solver,model2)
c1 = getobjval(model2)
@test c1 == Inf

solver.Lower_Prob = (X,k,p,opt,temp) -> [opt[1](X).lo,mid.(X),true,[]]
solver.Upper_Prob = (X,k,p,opt,temp) -> [opt[1](X).hi,mid.(X),false,[]]
model3 = BnBModel(X)
solveBnB!(solver,model3)
c2 = getobjval(model3)
@test c2 == Inf

solver.Lower_Prob = (X,k,p,opt,temp) -> [opt[1](X).lo,mid.(X),true,[]]
solver.Upper_Prob = (X,k,p,opt,temp) -> [opt[1](X).hi,mid.(X),true,[]]
function test_check(S,M,X1,X2)
    return (M.lbcnt == 2)
end
solver.Repeat_Chk = test_check
model5 = BnBModel(X)
solveBnB!(solver,model5)
c5 = getobjval(model5)
@test c5 == 1.5625

solver.Preprocess =  (feas_Pre,nsBox,UBDg,k_int,pos,opt,LBDn,UBDn,x,y) -> (false,nsBox)
solver.Lower_Prob = (X,k,p,opt,temp) -> [opt[1](X).lo,mid.(X),true,[]]
solver.Upper_Prob = (X,k,p,opt,temp) -> [((k==2) ? opt[1](X).hi : Inf),mid.(X),false,[]]
model4 = BnBModel(X)
solveBnB!(solver,model4)
c3 = getobjval(model4)
@test c3 == Inf
end
