"""
    Ipopt_UBD_eval_g!

Evaluates the constraint function g in place. Inputs:

* `x::Vector{Float64}`: Point to evaluate in X
* `g::Vector{Float64}`: Storage vector for resulting constraint eval.
* `opt::EAGO_Inner_NLP`: Option type containing problem information

No value returned. The function mutates g in place.
"""
function IPOPT_UBD_eval_g!(x::Vector{Float64},
                           g::Vector{Float64},
                           opt::EAGO_Inner_NLP)
    if opt.numConstr>0
        g[:] = opt.g(x)
    else
        g[:] = [-ones(x[1])]
    end
end

"""
    Ipopt_UBD_eval_grad_f!

Evaluates the gradient the objective f in place via Forward AD.
* `x::Vector{Float64}`: Point to evaluate in X
* `grad_f::Vector{Float64}`: Storage vector for resulting gradient of f.
* `opt`:Option type containing problem information

No value returned. The function mutates grad_f in place.
"""
function IPOPT_UBD_eval_grad_f!(x::Vector{Float64},grad_f::Vector{Float64},opt)
    ForwardDiff.gradient!(grad_f,opt.f,x)
end

"""
    Ipopt_UBD_eval_jac_g!


Evaluates the jacobian of the relaxation of the constraint function `g`. Inputs:
* `x::Vector{Float64}`: Point to evaluate in `X`
* `mode::Symbol`: Mode used for Jacobian/Hessian evaluation
* `rows::Vector{Int32}`: Row indices for nonzero Jacobian entries if structure
                         mode used
* `cols::Vector{Int32}`: Column Indices for nonzero Jacobian entries if
                         structure mode used
* `values::Vector{Float64}`: Storage for jacobian of convex relaxation of
                             constaint function g
* `opt_inner: inner_global`: Option type containing problem information

None. Mutates values in place.
"""
function IPOPT_UBD_eval_jac_g!(x::Vector{Float64},
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
            values[:] = transpose(ForwardDiff.jacobian(opt.g,x))
        else
            values[:] = zeros(x)
        end
    end
end

"""
    IPOPT_UBD_eval_h

Evaluates the hessian of the langrangian of the problem. Inputs:
* `x::Vector{Float64}`: Point to evaluate in X
* `mode::Symbol`: Mode used for Jacobian/Hessian evaluation
* `rows::Vector{Int32}`: Row indices for nonzero Jacobian entries if
                         structure mode used
* `cols::Vector{Int32}`: Column Indices for nonzero Jacobian entries if
                         structure mode used
* `obj_factor::Float64`: Scaling factor for objective
* `lambda::Vector{Float64}`: Multipliers for the constraints
* `values::Vector{Float64}`: Storage for langriangian hessian of convex
                            relaxation of the problem.
* `opts::EAGO_Inner_NLP`: Option type containing problem information

Returns nothing. Mutates values in place.
"""
function IPOPT_UBD_eval_h(x::Vector{Float64},
                          mode::Symbol,
                          rows::Vector{Int32},
                          cols::Vector{Int32},
                          obj_factor::Float64,
                          lambda::Vector{Float64},
                          values::Vector{Float64},
                          opt::EAGO_Inner_NLP)
    if mode == :Structure
      # Symmetric matrix, fill the lower left triangle only
      idx::Int64 = 1
      for row = 1:opt.numVar
        for col = 1:row
          rows[idx] = row
          cols[idx] = col
          idx += 1
        end
      end
    else
      if opt.numConstr>0
          c1::Array{Float64,2} = obj_factor*ForwardDiff.hessian(opt.f,x)+vector_hessian(opt.g, x, lambda)
      else
          c1 = obj_factor*ForwardDiff.hessian(opt.f,x)
      end
      val_arr::Array{Float64,1} = zeros(Float64,sum(1:opt.numVar))
      count::Int64 = 1
      offset::Int64 = 0
      for i=1:opt.numVar
          for j=(1+offset):(opt.numVar)
              val_arr[count] = c1[j,i]
              count += 1
          end
          offset += 1
      end
      values[:] = val_arr'
    end
end

"""
    Ipopt_UBD

Solves a upper bounding problem (local NLP) using Ipopt. Ipopt options and
domain reduction settings are controlled in the B&B algorithm using the
global_options type. Inputs:
* `X::Vector{Interval{Float64}}`: Node over which to solve the upper problem
* `k::Int64`: Number of Iterations The B&B Algorithm has taken
* `pos::Int64`: Depth of Node in B&B Tree
* `opt_inner`: Option type containing problem information
* `temp`: Info passed from lower to upper problem

Returns a tuple `(val,pnt,feas,X,Any[feas,val])` where
*`val::Float64`: Lower bound calculated
*`pnt::Array{Float64,1}`: An array of length equal to X that gives the
                          optimal solution of the lower bound problem.
* `feas::Bool`: Returns true if the problem is feasible and false if it is
                infeasible
* `Any[feas,val]`: Info for OBBT.
"""
function Ipopt_UBD(X,
                   k::Int64,
                   pos::Int64,
                   opts,
                   temp)

            # sets up problem
            x_L::Vector{Float64} = [X[i].lo for i=1:opts[1].numVar]
            x_U::Vector{Float64} = [X[i].hi for i=1:opts[1].numVar]

            prob = createProblem(opts[1].numVar, x_L, x_U, opts[1].numConstr, opts[1].gL, opts[1].gU,
                                 opts[1].numVar*opts[1].numConstr, Int64(opts[1].numVar*(opts[1].numVar+1)/2),
                                 x::Vector{Float64} -> opts[1].f(x),
                                 (x::Vector{Float64}, g::Vector{Float64}) -> opts[2].IPOPT_UBD_eval_g!(x,g),
                                 (x::Vector{Float64}, f::Vector{Float64}) -> opts[2].IPOPT_UBD_eval_grad_f!(x,f),
                                 (x::Vector{Float64}, mode::Symbol, rows::Vector{Int32}, cols::Vector{Int32}, values::Vector{Float64}) -> opts[2].IPOPT_UBD_eval_jac_g!(x, mode, rows, cols, values),
                                 (x::Vector{Float64}, mode::Symbol, rows::Vector{Int32}, cols::Vector{Int32}, obj_factor::Float64, lambda::Vector{Float64}, values::Vector{Float64}) -> opts[2].IPOPT_UBD_eval_h(x, mode, rows, cols, obj_factor, lambda, values))


            prob.x = mid.(X)

            if (opts[1].numConstr == 0)
                addOption(prob, "hessian_approximation", "limited-memory")
            end
            addOption(prob, "print_level", 0)

            if (opts[1].solver.UBD_full_depth < pos)
                addOption(prob, "tol", Inf)
            else
                addOption(prob, "tol", 1E-4)
            end
            # solve problem and unpacks variables
            TT = STDOUT
            redirect_stdout()
            status = solveProblem(prob)
            redirect_stdout(TT)
            pnt::Vector{Float64} = prob.x
            val::Float64 = prob.obj_val
            if (status == 0 || status == 1 || status == 6)
                feas::Bool = true
            elseif (status == 2)
                feas == false
            else
                error("Solver error code $status in Ipopt. Solution routine terminated.")
            end

            # output
            return val, pnt, feas, Any[feas,val]
end
