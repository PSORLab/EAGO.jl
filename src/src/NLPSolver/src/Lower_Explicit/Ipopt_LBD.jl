"""
    IPOPT_LBD_eval_f

Evaluates the convex relaxation of the objective, f. Inputs are:
* `x::Vector{Float64}`: Point to evaluate in X
* `X::Vector{Interval{Float64}}:` Node over which to solve the lower problem
* `opts`: Option type containing problem information
Returns the convex relaxation of the objective (::Float64).
"""
function IPOPT_LBD_eval_f(x::Vector{Float64},
                          X::Vector{Interval{Float64}},
                          opts::EAGO_Inner_NLP)
    f_SMC::HybridMC{opts.numVar,Interval{Float64},Float64} =  opts.f([HybridMC{opts.numVar,Interval{Float64},Float64}(SMCg{opts.numVar,Interval{Float64},Float64}(x[i],x[i],
                                                                          seed_g(Float64,i,opts.numVar),
                                                                          seed_g(Float64,i,opts.numVar),
                                                                          X[i],
                                                                          false)) for i=1:opts.numVar])
    return f_SMC.SMC.cv
end

"""
    IPOPT_LBD_eval_grad_f!

Evaluates the gradient of the convex relaxation of the objective f in place. Inputs:
* `x::Vector{Float64}`: Point to evaluate in X
* `X::Vector{Interval{Float64}}`: Node over which to solve the lower problem
* `f_grad:Vector{Float64}`: Storage vector for resulting gradient of f.
* `opts`: Option type containing problem information
No value returned. The function mutates `f_grad` in place.
"""
function IPOPT_LBD_eval_grad_f!(x::Vector{Float64},X::Vector{Interval{Float64}},
                                f_grad::Vector{Float64},opts::EAGO_Inner_NLP)

    f_SMC::HybridMC{opts.numVar,Interval{Float64},Float64} =  opts.f([HybridMC{opts.numVar,Interval{Float64},Float64}(SMCg{opts.numVar,Interval{Float64},Float64}(x[i],x[i],
                                                                          seed_g(Float64,i,opts.numVar),
                                                                          seed_g(Float64,i,opts.numVar),
                                                                          X[i],
                                                                          false)) for i=1:opts.numVar])

    f_grad[:] = f_SMC.SMC.cv_grad
end

"""
    IPOPT_LBD_eval_g!

Evaluates the convex relaxation of the constraint function g in place. Inputs:
* `x::Vector{Float64}`: Point to evaluate in X
* `X::Vector{Float64}`: Node over which to solve the lower problem
* `g::Vector{Float64}`: Storage vector for resulting constraint eval
* `opts`: Option type containing problem information
No value returned. The function mutates g in place.
"""
function IPOPT_LBD_eval_g!(x::Vector{Float64},
                           X::Vector{Interval{Float64}},
                           g::Vector{Float64},
                           opts::EAGO_Inner_NLP)
    if opts.numConstr>0
        x_SMC::Vector{HybridMC{opts.numVar,Interval{Float64},Float64}} = [HybridMC{opts.numVar,Interval{Float64},Float64}(SMCg{opts.numVar,Interval{Float64},Float64}(x[i],
                                                                              x[i],
                                                                              seed_g(Float64,i,opts.numVar),
                                                                              seed_g(Float64,i,opts.numVar),
                                                                              X[i],
                                                                              false)) for i=1:opts.numVar]

        g[:] = [opts.g(x_SMC)[i].SMC.cv for i=1:opts.numConstr]
    else
        g[:] = [-one(x[1])]
    end
end

"""
    IPOPT_LBD_eval_jac_g!

Evaluates the jacobian of the convex relaxation of the constraint function g. Inputs:
* `x::Vector{Float64}`: Point to evaluate in X
* `X::Vector{Interval{Float64}}`: Node over which to solve the lower problem
* `mode::Symbol`: Mode used for Jacobian/Hessian evaluation
* `rows::Vector{Int64}`: Row indices for nonzero Jacobian entries if structure
                         mode used
* `cols::Vector{Int64}`: Column Indices for nonzero Jacobian entries if
                         structure mode used
* `values::Array{Float64,2}` Storage for jacobian of convex relaxation of
                             constaint function g
* `opts`: Option type containing problem information
Returns nothing. Mutates values (::Array{Float64,2}) in place.
"""
function IPOPT_LBD_eval_jac_g!(Y,
                               x::Vector{Float64},
                               mode::Symbol,
                               rows::Vector{Int32},
                               cols::Vector{Int32},
                               values::Vector{Float64},
                               opt::EAGO_Inner_NLP,
                               cb::callback_storage)

    if mode == :Structure
       rows[:] = cb.row_temp_Ipopt_LBD
       cols[:] = cb.col_temp_Ipopt_LBD
    else
       if opt.numConstr>0

           x_SMC::Vector{HybridMC{opt.numVar,Interval{Float64},Float64}} = [HybridMC{opt.numVar,Interval{Float64},Float64}(SMCg{opt.numVar,Interval{Float64},Float64}(x[i],
                                                                             x[i],
                                                                             seed_g(Float64,i,opt.numVar),
                                                                             seed_g(Float64,i,opt.numVar),
                                                                             Y[i],
                                                                             false)) for i=1:opt.numVar]
           g_SMC::Vector{HybridMC{opt.numVar,Interval{Float64},Float64}} = opt.g(x_SMC)
           g_jac::Array{Float64,2} = [g_SMC[i].SMC.cv_grad[j] for j=1:opt.numVar, i=1:opt.numConstr]
           values[:] = g_jac
       else
           values[:] = zeros(x)
     end
    end
end

"""
    Ipopt_LBD

Solves a lower bounding problem (NLP) constructed by relaxing the subproblem
problem using the twice-differentiable McCormick operators presented by Khan2017.
The relaxed problem is then solved using Ipopt. Ipopt options and domain
reduction settings are controlled in the B&B algorithm using the global_options
type. Inputs:
* `X::Vector{Interval{Float64}}`: Node over which to solve the lower problem
* `k::Int64`: Number of Iterations The B&B Algorithm has taken
* `pos::Int64`: Depth of Node in B&B Tree
* `opt`: Option type containing problem information
* `UBD::Float64`: Global upper bound for B&B algorithm
Returns a tuple `(val,pnt,feas,X,[])` where
* `val::Float64` Lower bound calculated
* `pnt::Array{Float64,1}`: An array of length equal to X that gives the
                           optimal solution of the lower bound problem.
* `feas::Bool` Returns true if the problem is feasible and false if it is
               infeasible
* `temp`: Information useful for OBBT, Any[mult_lo,mult_hi,val]
"""
function Ipopt_LBD(X::Vector{Interval{Float64}},k::Int64,pos::Int64,opts::Any,UBD::Float64)

    opts[1].solver.SubGradRefine && set_hybrid_box!(SVector{opts[1].numVar,Interval{Float64}}(X),SVector{opts[1].numVar,Float64}(mid.(X)),true)

    # sets up problem
    x_L::Vector{Float64} = [X[i].lo for i=1:opts[1].numVar]
    x_U::Vector{Float64} = [X[i].hi for i=1:opts[1].numVar]

    feqn = x::Vector{Float64} -> opts[2].IPOPT_LBD_eval_f(x,X)
    geqn = (x::Vector{Float64}, g::Vector{Float64}) -> opts[2].IPOPT_LBD_eval_g!(x,X,g)
    fgreqn = (x::Vector{Float64}, f::Vector{Float64}) -> opts[2].IPOPT_LBD_eval_grad_f!(x,X,f)
    gjaceqn = (x::Vector{Float64}, mode::Symbol, rows::Vector{Int32}, cols::Vector{Int32}, values::Vector{Float64}) -> opts[2].IPOPT_LBD_eval_jac_g!(X, x, mode, rows, cols, values)

    prob = createProblem(opts[1].numVar, x_L, x_U, opts[1].numConstr, opts[1].gL, opts[1].gU, opts[1].numVar*opts[1].numConstr, Int64(opts[1].numVar*(opts[1].numVar+1)/2),
                     feqn, geqn, fgreqn, gjaceqn)
    prob.x = mid.(X)

    addOption(prob, "print_level", 0)
    addOption(prob, "hessian_approximation", "limited-memory")

    # solve problem and unpacks variables
    status::Int64 = solveProblem(prob)
    pnt::Vector{Float64} = prob.x
    val::Float64 = prob.obj_val
    if (status == 0 || status == 1 || status == 6)
        feas = true
    elseif (status == 2)
        feas = false
    else
        error("Solver error code $status in Ipopt. Solution routine terminated.")
    end

    mult_lo::Vector{Float64} = prob.mult_x_L
    mult_hi::Vector{Float64} = prob.mult_x_U
    temp = Any[mult_lo,mult_hi,val]

    # Formats output and returns appropriate values
    return val, pnt, feas, temp
end
