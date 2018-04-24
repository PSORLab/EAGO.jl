#=
function aBB_IPOPT_LBD_eval_h(x::Vector{Float64},
                              X::Vector{Interval{Float64}},
                              mode::Symbol,
                              rows::Array{Int32,1},
                              cols::Array{Int32,1},
                              obj_factor::Float64,
                              lambda::Vector{Float64},
                              values::Array{Float64,1},
                              opts::EAGO_Inner_NLP,
                              cb::callback_storage,
                              finput::Function)
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
         c1::Array{Float64,2} = obj_factor*ForwardDiff.hessian(y::Vector{Float64} -> finput(y),x)+vector_hessian(y::Vector{Float64} -> cb.IPOPT_LBD_eval_g(y,X), x, lambda)
      else
         c1 = obj_factor*ForwardDiff.hessian(y::Vector{Float64} -> finput(y) ,x)
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

function AlphaBB_LBD(X::Vector{Interval{Float64}},k::Int64,pos::Int64,opt,UBD::Float64)
    ReverseDiff.hessian!(opt[2].alphaBB_Result,opt[2].alphaBB_Tape,X)
    alpha::Vector{Float64} = zeros(Float64,opt[1].numVar)
    for i=1:opt[1].numVar
        hcalc::Float64 = opt[2].alphaBB_Result[i,i].lo
        for j=1:opt[1].numVar
            if (i != j)
                hcalc += -max(abs(opt[2].alphaBB_Result[i,j].lo),
                              abs(opt[2].alphaBB_Result[i,j].hi))*(X[j].hi-X[j].lo)/(X[i].hi-X[i].lo)
            end
        end
        alpha[i] = max(0.0,-0.5*hcalc)
    end

    # feasibility-based bound tightening & Domain Reduction Presolves
    feas::Bool = true
    composite_DR_pre!(X,feas,k,pos,opt,UBD)

    # Solve the relaxed optimization problem
    if (feas == true)
        # sets up problem
        x_L::Vector{Float64} = [X[i].lo for i=1:opt[1].numVar]
        x_U::Vector{Float64} = [X[i].hi for i=1:opt[1].numVar]
        prob = createProblem(opt[1].numVar, x_L, x_U, opt[1].numConstr, opt[1].gL, opt[1].gU, opt[1].numVar*opt[1].numConstr, Int64(opt[1].numVar*(opt[1].numVar+1)/2),
                             x::Vector{Float64} -> opt[1].f(x) - sum([alpha[i]*(x-X[i].lo)*(X[i].hi-x) for i=1:opt[1].numVar]),
                            (x::Vector{Float64}, g::Vector{Float64}) -> opt[2].IPOPT_LBD_eval_g!(x,X,g),
                            (x::Vector{Float64}, f::Vector{Float64}) -> ForwardDiff.gradient!(f,x -> opt[1].f(x) - sum([alpha[i]*(x-X[i].lo)*(X[i].hi-x) for i=1:opt[1].numVar]),x),
                            (x::Vector{Float64}, mode::Symbol, rows::Array{Int32,1}, cols::Array{Int32,1}, values::Vector{Float64}) -> opt[2].IPOPT_LBD_eval_jac_g!(x, X, mode, rows, cols, values),
                            (x::Vector{Float64}, mode::Symbol, rows::Array{Int32,1}, cols::Array{Int32,1}, obj_factor::Float64, lambda::Vector{Float64}, values::Vector{Float64})
                                -> opt[2].aBB_IPOPT_LBD_eval_h(x,X,mode,rows,cols,obj_factor,lambda,values,opt[1],opt[2],
                                   x::Vector{Float64} -> opt[1].f(x) - sum([alpha[i]*(x-X[i].lo)*(X[i].hi-x) for i=1:opt[1].numVar])))
        prob.x = mid.(X)
        addOption(prob, "hessian_approximation", "limited-memory")
        addOption(prob, "tol", 1E-5)
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
    end

    mult_lo::Vector{Float64} = [0.0 for i=1:opt[1].numVar]
    mult_hi::Vector{Float64} = [0.0 for i=1:opt[1].numVar]
    for i=1:opt[1].numVar
        if abs(mult[i])>opt[1].solver.dual_tol
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

    # Formats Runs Post Contractor
    composite_DR_post!(X,k,pos,opt[1],val,UBD,mult_lo,mult_hi)

    # Formats output and returns appropriate values
    return val, pnt, feas, X, []
end
=#
