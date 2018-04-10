workspace()

using IntervalArithmetic
using EAGO

X = [Interval(1,5)]
model = BnBModel(X)
solver = BnBSolver()
set_to_default!(solver)
solver.opt = [x -> x[1]^2]
solver.Lower_Prob = (X,k,p,opt,temp) -> [opt[1](X).lo,mid.(X),true,[]]
solver.Upper_Prob = (X,k,p,opt,temp) -> [opt[1](X).hi,mid.(X),true,[]]
a = solver.Lower_Prob(X,1,1,solver.opt,[])
b = solver.Upper_Prob(X,1,1,solver.opt,[])
solveBnB!(solver,model)

solver.iter_lim = true
solver.max_iter = 10
solver.Lower_Prob = (X,k,p,opt,temp) -> [opt[1](X).lo,mid.(X),false,[]]
model2 = BnBModel(X)
solveBnB!(solver,model2)
c1 = getobjval(model2)

solver.Lower_Prob = (X,k,p,opt,temp) -> [opt[1](X).lo,mid.(X),true,[]]
solver.Upper_Prob = (X,k,p,opt,temp) -> [opt[1](X).hi,mid.(X),false,[]]
model3 = BnBModel(X)
solveBnB!(solver,model3)
c2 = getobjval(model3)

solver.Lower_Prob = (X,k,p,opt,temp) -> [opt[1](X).lo,mid.(X),true,[]]
solver.Upper_Prob = (X,k,p,opt,temp) -> [((k==2) ? opt[1](X).hi : Inf),mid.(X),false,[]]
model4 = BnBModel(X)
solveBnB!(solver,model4)
c3 = getobjval(model4)
println("c3: $c3")
#@test c3 == Inf

solver.Lower_Prob = (X,k,p,opt,temp) -> [opt[1](X).lo,mid.(X),true,[]]
solver.Upper_Prob = (X,k,p,opt,temp) -> [opt[1](X).hi,mid.(X),true,[]]
solver.Repeat_Chk = (S,M,X1,X2) -> (M.lbcnt == 3)
model5 = BnBModel(X)
solveBnB!(solver,model5)
c5 = getobjval(model5)
println("c5: $c5")
#@test c4 == Inf
