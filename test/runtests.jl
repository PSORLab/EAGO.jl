#!/usr/bin/env julia

using EAGO, JuMP, MathOptInterface #, Ipopt
const MOI = MathOptInterface

println("BEGIN TESTING BRANCH AND BOUND LIBRARY...")
#include("branch_bound.jl")
println("TESTING BRANCH AND BOUND LIBRARY COMPLETE.")

println("BEGIN TESTING NODE ACCESS ROUTINES...")
include("node_access.jl")                                # WORKING
println("TESTING NODE ACCESS ROUTINES COMPLETE.")

println("BEGIN TESTING PARAMETRIC INTERVAL LIBRARY...")
#include("ParametricInterval/runtests.jl")
println("TESTING PARAMETRIC INTERVAL LIBRARY COMPLETE.")

println("BEGIN TESTING MCCORMICK LIBRARY...")
#nclude("McCormick/mccormick.jl")                         # WORKING EXCEPT REVERSE
println("TESTING MCCORMICK LIBRARY COMPLETE.")

println("BEGIN DOMAIN REDUCTION LIBRARY...")
include("domain_reduction.jl")                             # WORKING
println("TESTING DOMAIN REDUCTION COMPLETE.")

println("BEGIN RELAXATION ROUTINES...")
include("relaxations.jl")                                 # WORKING
println("TESTING RELAXATION ROUTINES COMPLETE.")

println("BEGIN TESTING OPTIMIZER")
include("Optimizer/optimizer.jl")                        # WORKING
println("END TESTING OPTIMIZER")

println("BEGIN TESTING SCRIPT BRIDGE")
include("script_optimizer.jl")                            # WORKING
println("END TESTING SCRIPT BRIDGE")

println("BEGIN TESTING SIP SOLVER")
#include("semiinfinite.jl")
println("END TESTING SIP SOLVER")


#=
low_eval =  EAGO.ImplicitLowerEvaluator{1}()
new_node = EAGO.NodeBB(Float64[-0.78,6.0],Float64[-0.40,9.0],-Inf,Inf,0,-1,false)
ytest = [-0.59, 7.5]
ytest2 = [-0.59,6.5]
EAGO.set_current_node!(low_eval, new_node)
test1 = low_eval.current_node.lower_variable_bounds
test2 = low_eval.current_node.upper_variable_bounds

f(x,p) = x[1]
g(x,p) = [x[1]]
function h(out,x,p)
    t1 = x[1]^2
    t2 = x[1]*p[1]
    t3 = 4.0
    t4 = t1 + t2
    t5 = t4 + t3
    out[1] = t5
end

function hj(out,x,p)
    out[1] = 2*x[1]+p[1]
end
np = 1; nx = 1; ng = 1
sparse_pattern = Tuple{Int64,Int64}[(1,1)]

EAGO.build_evaluator!(low_eval, h, np, nx)
EAGO.build_evaluator!(low_eval, h, np, nx, obj = f)
EAGO.build_evaluator!(low_eval, h, np, nx, obj = f, constr = g, ng = 1)
EAGO.build_evaluator!(low_eval, h, np, nx, obj = f, constr = g, ng = 1)
EAGO.build_evaluator!(low_eval, h, np, nx, obj = f, constr = g, ng = 1, state_jac = hj)
EAGO.build_evaluator!(low_eval, h, np, nx, obj = f, constr = g, ng = 1, state_jac = hj, user_sparse = sparse_pattern)

test3 = low_eval.state_relax
test4 = low_eval.cnstr_relax
test5 = low_eval.state_ref_relaxation
test6 = isnan(low_eval.last_p[1])
test7 = low_eval.ref_p

# initial implicit eval
EAGO.relax_implicit!(low_eval, ytest)

test8 = low_eval.obj_eval
test9 = low_eval.cnstr_eval
test10 = low_eval.ref_p
test11 = isnan(low_eval.last_p[1])
test12 = low_eval.P
test13 = low_eval.X

test14 = low_eval.state_relax
test15 = low_eval.state_ref_relaxation

# second implicit eval
EAGO.relax_implicit!(low_eval, ytest2)
test16 = low_eval.obj_eval
test17 = low_eval.cnstr_eval
test18 = low_eval.ref_p[1]
test19 = isnan(low_eval.last_p[1])
test20 = low_eval.state_relax
test21 = low_eval.state_ref_relaxation

# test objective eval
EAGO.relax_objective!(low_eval)
test31 = low_eval.obj_relax
test32 = low_eval.obj_eval

# test constraint eval
EAGO.relax_constraints!(low_eval)
test41 = low_eval.cnstr_relax
test42 = low_eval.cnstr_eval

test51 = MOI.eval_objective(low_eval, ytest)

gval = zeros(1)
MOI.eval_constraint(low_eval, gval, ytest)
test61 = gval[1]

df = zeros(1)
MOI.eval_objective_gradient(low_eval, df, ytest)
test71 = df[1]

test81 = MOI.jacobian_structure(low_eval)

jac_val = zeros(1,1)
MOI.eval_constraint_jacobian(low_eval, jac_val, ytest)
test91 = jac_val[1]

w = 0.5
jac_prod_val = zeros(1)
MOI.eval_constraint_jacobian_product(low_eval, jac_prod_val, ytest, w)
test10_1 = jac_prod_val[1]

jact_prod_val = zeros(1)
MOI.eval_constraint_jacobian_transpose_product(low_eval, jact_prod_val, ytest, w)
test10_2 = jact_prod_val[1]

features = MOI.features_available(low_eval)
test11_1 = features[1]
test11_2 = features[2]
=#
