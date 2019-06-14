using JuMP, Ipopt

f(x) = x[1]^2/3.0 + x[2]^2 + x[1]/2.0
g(x,y) = (1.0 - (x[1]^2)*(y[1]^2))^2 - x[1]*y[1]^2 - x[2]^2 + x[2]

m = Model(with_optimizer(Ipopt.Optimizer))
@variable(m, x[1:5])
gi = :g1
eval(:(($gi)(x...) = g(x,1.5)))
eval(:(register(m, :($gi), 5, $gi, autodiff=true)))
func_call = Expr(:call)
args = []
push!(args, gi)
for i in 1:nx
    push!(args, JuMP.VariableRef(m,i))
end
ineq_call = Expr(:call, [func_call; :(<=), -eps_g])
add_NL_constraint(m, ineq_call)

#add_NL_constraint(m, :(g1(x...) <= -0.7))
#user_ops = m.nlp_data.user_operators
#add_NL_constraint(m, :(g1(x...) <= -0.9))
#optimize!(m)


#=
optimizer = backend(m)#.optimizer
evaluator = optimizer.optimizer.model.optimizer.nlp_data.evaluator
mJuMP = evaluator.m
bmJuMP = backend(mJuMP)
=#
