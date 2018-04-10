"""
    snopt_callback_LBD

Callback function used by SNOPT in a lower bounding problem (NLP) constructed by
relaxing the subproblem using the differentiable McCormick operators presented
by Khan2017. Inputs:
* `y::Vector{Float64}`  - Point to evaluate in X
* `X::Vector{Interval{Float64}}` - Node corresponding to lower problem
* `opts` - Option type containing problem information
Returns a tuple tuple (f_cv, c_cv, dfdx_cv, dcdx_cv, fail) where
* `f_cv::Float64`: Convex relaxation of the objective function
* `c_cv::Array{Float64,1}`: Convex relaxation of the constraint function
* `dfdx_cv::Array{Float64,1}`: Gradient of the convex relaxation of the
                               objective function
* `dcdx_cv::Array{Float64,2}`: Jacobian of the convex relaxation of the
                               constraint function
* `fail::Bool`: Flag used by SNOPT solver to indicate if SNOPT solver failed.
"""
function snopt_callback_LBD(y::Vector{Float64},
                            X::Vector{Interval{Float64}},
                            opt)

    # Create SMCg relaxtion of y on X
    x_SMC::Vector{SMCg{opt.numVar,Interval{Float64},Float64}} = [SMCg{opt.numVar,Interval{Float64},Float64}(y[i],
                                                                        y[i],
                                                                        seed_g(Float64,i,opt.numVar),
                                                                        seed_g(Float64,i,opt.numVar),
                                                                        X[i],
                                                                        false,
                                                                        SVector{opt.numVar,Interval{Float64}}(X),
                                                                        SVector{opt.numVar,Float64}(y)) for i=1:opt.numVar]

    # Evaluates the relaxation of objective and constraints
    # relaxation of function
    f_val::SMCg{opt.numVar,Interval{Float64},Float64} = opt.f(x_SMC)
    #println("x_SMC: ",x_SMC)
    #println("f_val",f_val)


    if opt.numConstr>0
        c::Vector{SMCg{opt.numVar,Interval{Float64},Float64}} = opt.g(x_SMC)
        dcdx::SparseMatrixCSC{Float64,Int64} = spzeros(length(opt.gL_loc)+length(opt.gU_loc),opt.numVar)
    else
        dcdx = spzeros(1,opt.numVar)
    end

    # forms coefficients of linear relaxations
    if opt.numConstr>0
        cx_ind1::Int64 = 1
        for i in opt.gL_loc
            for j=1:opt.numVar
                if (c[i].cv_grad[j] != 0.0)
                    dcdx[cx_ind1,j] = c[i].cv_grad[j]
                end
            end
            cx_ind1 += 1
        end
        for i in opt.gU_loc
            for j=1:opt.numVar
                if (c[i].cc_grad[j] != 0.0)
                    dcdx[cx_ind1,j] = -c[i].cc_grad[j]
                end
            end
            cx_ind1 += 1
        end
    end

    # forms rhs of linear relaxations
    if opt.numConstr>0
        rhs::Vector{Float64} = zeros(Float64,length(opt.gL_loc)+length(opt.gU_loc))
    else
        rhs = zeros(Float64,1)
    end

    if opt.numConstr>0
        cx_ind2::Int64 = 1
        for i in opt.gU_loc
            rhs[cx_ind2] = c[i].cv - opt.gU[i]
            cx_ind2 += 1
        end
        for i in opt.gL_loc
            rhs[cx_ind2] = opt.gL[i]-c[i].cc
            cx_ind2 += 1
        end
    end

    fail::Bool = false
    return f_val.cv, rhs, f_val.cv_grad, dcdx, fail
end

function snopt_callback_LBD(y::Vector{Float64},
                            X::Vector{MCInterval{Float64}},
                            opt)

    # Create SMCg relaxtion of y on X
    x_SMC::Vector{SMCg{opt.numVar,MCInterval{Float64},Float64}} = [SMCg{opt.numVar,MCInterval{Float64},Float64}(y[i],
                                                                        y[i],
                                                                        seed_g(Float64,i,opt.numVar),
                                                                        seed_g(Float64,i,opt.numVar),
                                                                        X[i],
                                                                        false,
                                                                        SVector{opt.numVar,MCInterval{Float64}}(X),
                                                                        SVector{opt.numVar,Float64}(y)) for i=1:opt.numVar]

    # Evaluates the relaxation of objective and constraints
    # relaxation of function
    f_val::SMCg{opt.numVar,MCInterval{Float64},Float64} = opt.f(x_SMC)
    #println("x_SMC: ",x_SMC)
    #println("f_val",f_val)


    if opt.numConstr>0
        c::Vector{SMCg{opt.numVar,MCInterval{Float64},Float64}} = opt.g(x_SMC)
        dcdx::SparseMatrixCSC{Float64,Int64} = spzeros(length(opt.gL_loc)+length(opt.gU_loc),opt.numVar)
    else
        dcdx = spzeros(1,opt.numVar)
    end

    # forms coefficients of linear relaxations
    if opt.numConstr>0
        cx_ind1::Int64 = 1
        for i in opt.gL_loc
            for j=1:opt.numVar
                if (c[i].cv_grad[j] != 0.0)
                    dcdx[cx_ind1,j] = c[i].cv_grad[j]
                end
            end
            cx_ind1 += 1
        end
        for i in opt.gU_loc
            for j=1:opt.numVar
                if (c[i].cc_grad[j] != 0.0)
                    dcdx[cx_ind1,j] = -c[i].cc_grad[j]
                end
            end
            cx_ind1 += 1
        end
    end

    # forms rhs of linear relaxations
    if opt.numConstr>0
        rhs::Vector{Float64} = zeros(Float64,length(opt.gL_loc)+length(opt.gU_loc))
    else
        rhs = zeros(Float64,1)
    end

    if opt.numConstr>0
        cx_ind2::Int64 = 1
        for i in opt.gU_loc
            rhs[cx_ind2] = c[i].cv - opt.gU[i]
            cx_ind2 += 1
        end
        for i in opt.gL_loc
            rhs[cx_ind2] = opt.gL[i]-c[i].cc
            cx_ind2 += 1
        end
    end

    fail::Bool = false
    return f_val.cv, rhs, f_val.cv_grad, dcdx, fail
end


"""
    SNOPT_LBD

Solves a lower bounding problem (NLP) constructed by relaxing the subproblem
problem using the differentiable McCormick operators presented by Khan2017. The
relaxed problem is then solved using SNOPT. SNOPT options and domain reduction
settings are controlled in the B&B algorithm using the global_options type. Inputs:
* `X::Vector{Interval{Float64}}`: Node over which to solve the lower problem
* `k::Int64`: Number of Iterations The B&B Algorithm has taken
* `pos::Int64`: Depth of Node in B&B Tree
* `opt::Any`: Option type containing problem information (a solver type)
* `UBD::Float64` - Global upper bound for B&B algorithm
Returns a tuple `(val,pnt,feas,X,Any[mult_lo,mult_hi,val])` where
* `val::Float64`: Lower bound calculated
* `pnt::Array{Float64,1}`: An array of length equal to X that gives the
                           optimal solution of the lower bound problem.
* `feas::Bool`: Returns true if the problem is feasible and false if it is
                infeasible
* `Any[mult_lo,mult_hi,val]`: Information for use in optimality-based bound tightening.
"""
function SNOPT_LBD(X::Vector{Interval{Float64}},
                   k::Int64,
                   pos::Int64,
                   opt,
                   UBD::Float64)
        x_Li::Vector{Float64} = [X[i].lo for i=1:opt[1].numVar]
        x_Ui::Vector{Float64} = [X[i].hi for i=1:opt[1].numVar]
        x0i::Vector{Float64} = (x_Li + x_Ui)/2.0
        options = Dict{String, Any}()
        options["Function precision"] = 1.0E-12
        options["Derivative option"] = 1
        options["Print file"] = 0
        options["Major print level"] = 0
        options["Minor print level"] = 0
        options["Major optimality tolerance"] = 1e-6
        options["Verify level"] = 0

        #TT = STDOUT
        #redirect_stdout()
        pnt::Vector{Float64}, val::Float64, status::Int64, mult::Vector{Float64} = Snopt.snopt(y::Vector{Float64} -> opt[2].fg_SNOPT_LBD(y,X), x0i, x_Li, x_Ui, options)
        #redirect_stdout(TT)
        #tups = snopt_callback(pnt,X)
        #val = tups[1]
        #println("snopt_callback(x,X)",snopt_callback(x,X))
        if (status == 1 || status == 2)
            feas = true
        elseif (status == 11 || status == 12 || status == 13 || status == 14)
            feas = false
        else
            error("Solver error code $status in Snopt. Solution routine terminated.")
        end

        # Optimality-based bound tightening on variables only
        if (opt[1].solver.variable_depth >= pos)
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
        else
            mult_lo = Float64[]
            mult_hi = Float64[]
        end
        temp = Any[mult_lo,mult_hi,val]
        # Formats output and returns appropriate values
        return val, pnt, feas, temp
end

function SNOPT_LBD(X::Vector{MCInterval{Float64}},
                   k::Int64,
                   pos::Int64,
                   opt,
                   UBD::Float64)
        x_Li::Vector{Float64} = [X[i].lo for i=1:opt[1].numVar]
        x_Ui::Vector{Float64} = [X[i].hi for i=1:opt[1].numVar]
        x0i::Vector{Float64} = (x_Li + x_Ui)/2.0
        options = Dict{String, Any}()
        options["Function precision"] = 1.0E-12
        options["Derivative option"] = 1
        options["Print file"] = 0
        options["Major print level"] = 0
        options["Minor print level"] = 0
        options["Major optimality tolerance"] = 1e-6
        options["Verify level"] = 0

        #TT = STDOUT
        #redirect_stdout()
        pnt::Vector{Float64}, val::Float64, status::Int64, mult::Vector{Float64} = Snopt.snopt(y::Vector{Float64} -> opt[2].fg_SNOPT_LBD(y,X), x0i, x_Li, x_Ui, options)
        #redirect_stdout(TT)
        #tups = snopt_callback(pnt,X)
        #val = tups[1]
        #println("snopt_callback(x,X)",snopt_callback(x,X))
        if (status == 1 || status == 2)
            feas = true
        elseif (status == 11 || status == 12 || status == 13 || status == 14)
            feas = false
        else
            error("Solver error code $status in Snopt. Solution routine terminated.")
        end

        # Optimality-based bound tightening on variables only
        if (opt[1].solver.variable_depth >= pos)
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
        else
            mult_lo = Float64[]
            mult_hi = Float64[]
        end
        temp = Any[mult_lo,mult_hi,val]
        # Formats output and returns appropriate values
        return val, pnt, feas, temp
end
