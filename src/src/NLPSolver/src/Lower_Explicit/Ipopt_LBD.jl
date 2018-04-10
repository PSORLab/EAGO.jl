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
    f_SMC::SMCg{opts.numVar,Float64} =  opts.f([SMCg{opts.numVar,Float64}(x[i],x[i],
                                                                          seed_g(Float64,i,opts.numVar),
                                                                          seed_g(Float64,i,opts.numVar),
                                                                          X[i],
                                                                          false,
                                                                          X,
                                                                          x) for i=1:opts.numVar])
    return f_SMC.cv
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
    f_SMC::SMCg{opts.numVar,Float64} =  opts.f([SMCg{opts.numVar,Float64}(x[i],x[i],
                                                                          seed_g(Float64,i,opts.numVar),
                                                                          seed_g(Float64,i,opts.numVar),
                                                                          X[i],
                                                                          false,
                                                                          X,
                                                                          x) for i=1:opts.numVar])
    f_grad[:] = f_SMC.cv_grad
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
        x_SMC::Vector{SMCg{opts.numVar,Float64}} = [SMCg{opts.numVar,Float64}(x[i],
                                                                              x[i],
                                                                              seed_g(Float64,i,opts.numVar),
                                                                              seed_g(Float64,i,opts.numVar),
                                                                              X[i],
                                                                              false,
                                                                              X,
                                                                              x) for i=1:opts.numVar]
        g[:] = [opts.g(x_SMC)[i].cv for i=1:opts.numConstr]
    else
        g[:] = -ones(x[1])
    end
end

"""
    IPOPT_LBD_eval_g

Evaluates the convex relaxation of the constraint function g.
* `x::Vector{Float64}`: Point to evaluate in X
* `X::Vector{Interval{Float64}}: Node over which to solve the lower problem
* `opts`: Option type containing problem information
Returns the convex relaxation of the constraint function g (::Vector{Float64}).
"""
function IPOPT_LBD_eval_g(x::Vector{Float64},
                          X::Vector{Interval{Float64}},
                          opts::EAGO_Inner_NLP)
    if opts.numConstr>0
        x_SMC::Vector{SMCg{opts.numVar,Float64}} = [SMCg{opts.numVar,Float64}(x[i],
                                                                              x[i],
                                                                              seed_g(Float64,i,opts.numVar),
                                                                              seed_g(Float64,i,opts.numVar),
                                                                              X[i],
                                                                              false,
                                                                              X,
                                                                              x) for i=1:opts.numVar]
        return [opts.g(x_SMC)[i].cv for i=1:opts.numConstr]
    else
        return [-ones(x[1])]
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
function IPOPT_LBD_eval_jac_g!(x::Vector{Float64},
                               X::Vector{Interval{Float64}},
                               mode::Symbol,
                               rows::Array{Int32,1},
                               cols::Array{Int32,1},
                               values::Array{Float64,1},
                               opts::EAGO_Inner_NLP,
                               cb::callback_storage)

    if mode == :Structure
        rows[:] = cb.col_temp_Ipopt_LBD
        cols[:] = cb.row_temp_Ipopt_LBD
    else
        x_SMC::Vector{SMCg{opts.numVar,Float64}} = [SMCg{opts.numVar,Float64}(x[i],
                                                                              x[i],
                                                                              seed_g(Float64,i,opts.numVar),
                                                                              seed_g(Float64,i,opts.numVar),
                                                                              X[i],
                                                                              false,
                                                                              X,
                                                                              x) for i=1:opts.numVar]
        g_SMC::Vector{SMCg{opts.numVar,Float64}} = opts.g(x_SMC)
        g_jac::Array{Float64,2} = [g_SMC[i].cv_grad[j] for i=1:opts.numConstr, j=1:opts.numVar]
        values[:] = transpose(g_jac)
    end
end

"""
    IPOPT_LBD_eval_h

Evaluates the hessian of the langrangian of convex relaxation of the problem. Inputs are:
* `x::Vector{Float64}`: Point to evaluate in X
* `X::Vector{Interval{Float64}}`: Node over which to solve the lower problem
* `mode::Symbol`: Mode used for Jacobian/Hessian evaluation
* `rows::Vector{Int64}`: Row indices for nonzero Jacobian entries if
                         structure mode used
* `cols::Vector{Int64}`: Column Indices for nonzero Jacobian entries if
                         structure mode used
* `obj_factor::Float64`: Scaling factor for objective
* `lambda::Vector{Float64}`: Multipliers for the constraints
* `values::Array{Float64,2}`: Storage for langriangian hessian of convex
                              relaxation of the problem.
* `opts`: Option type containing problem information
Returns nothing. Mutates values (::Array{Float64,2}) in place.
"""
function IPOPT_LBD_eval_h(x::Vector{Float64},
                          X::Vector{Interval{Float64}},
                          mode::Symbol,
                          rows::Array{Int32,1},
                          cols::Array{Int32,1},
                          obj_factor::Float64,
                          lambda::Vector{Float64},
                          values::Array{Float64,1},
                          opts::EAGO_Inner_NLP,
                          cb::callback_storage)
    if mode == :Structure
      # Symmetric matrix, fill the lower left triangle only
      idx::Int64 = 1
      for row = 1:opts.numVar
        for col = 1:row
          rows[idx] = cb.row_temp_Ipopt_LBD[idx]
          cols[idx] = cb.col_temp_Ipopt_LBD[idx]
          idx += 1
        end
      end
    else
      if opts.numConstr>0
         c1::Array{Float64,2} = obj_factor*ForwardDiff.hessian(y::Vector{Float64} -> cb.IPOPT_LBD_eval_f(y,X),x)+vector_hessian(y::Vector{Float64} -> cb.IPOPT_LBD_eval_g(y,X), x, lambda)
      else
         c1 = obj_factor*ForwardDiff.hessian(y::Vector{Float64} -> cb.IPOPT_LBD_eval_f(y,X),x)
     end
     val_arr::Vector{Float64} = zeros(Float64,sum(1:opts.numVar))
     count::Int64 = 1
     offset::Int64 = 0
     for i=1:opts.numVar
        for j=(1+offset):(opts.numVar)
            val_arr[count] = c1[j,i]
            count += 1
        end
        offset += 1
     end
     values[:] = val_arr'
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

    # sets up problem
    x_L::Vector{Float64} = [X[i].lo for i=1:opts[1].numVar]
    x_U::Vector{Float64} = [X[i].hi for i=1:opts[1].numVar]
    prob = createProblem(opts[1].numVar, x_L, x_U, opts[1].numConstr, opts[1].gL, opts[1].gU, opts[1].numVar*opts[1].numConstr, Int64(opts[1].numVar*(opts[1].numVar+1)/2),
                     x::Vector{Float64} -> opts[2].IPOPT_LBD_eval_f(x,X),
                     (x::Vector{Float64}, g::Vector{Float64}) -> opts[2].IPOPT_LBD_eval_g!(x,X,g),
                     (x::Vector{Float64}, f::Vector{Float64}) -> opts[2].IPOPT_LBD_eval_grad_f!(x,X,f),
                     (x::Vector{Float64}, mode::Symbol, rows::Array{Int32,1}, cols::Array{Int32,1}, values::Vector{Float64}) -> opts[2].IPOPT_LBD_eval_jac_g!(x, X, mode, rows, cols, values),
                     (x::Vector{Float64}, mode::Symbol, rows::Array{Int32,1}, cols::Array{Int32,1}, obj_factor::Float64, lambda::Vector{Float64}, values::Vector{Float64}) -> opts[2].IPOPT_LBD_eval_h(x, X, mode, rows, cols, obj_factor, lambda, values))
    prob.x = mid.(X)

    addOption(prob, "hessian_approximation", "limited-memory")
    addOption(prob, "tol", 1E-6)
    addOption(prob, "print_level", 0)

    # solve problem and unpacks variables
    TT = STDOUT
    #redirect_stdout()
    status::Int64 = solveProblem(prob)
    #redirect_stdout(TT)
    pnt::Vector{Float64} = prob.x
    val::Float64 = prob.obj_val
    if (status == 0 || status == 1 || status == 6)
        feas = true
    elseif (status == 2)
        feas == false
    else
        error("Solver error code $status in Ipopt. Solution routine terminated.")
    end

    mult_lo::Vector{Float64} = [0.0 for i=1:opts[1].numVar]
    mult_hi::Vector{Float64} = [0.0 for i=1:opts[1].numVar]
    for i=1:opts[1].numVar
        if abs(mult[i])>opts[1].solver.dual_tol
            if (pnt[i]-X[i].lo)==(X[i].hi-pnt[i])
                mult_lo[i] = mult[i]
                mult_hi[i] = mult[i]
            elseif (pnt[i]-X[i].lo)<(X[i].hi-pnt[i])
                mult_lo[i] = mult[i]
            else
                mult_hi[i] = mult[i]
            end
        end
    end
    temp = Any[mult_lo,mult_hi,val]

    # Formats output and returns appropriate values
    return val, pnt, feas, temp
end
