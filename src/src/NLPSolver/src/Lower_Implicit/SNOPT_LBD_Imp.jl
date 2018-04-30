"""
    snopt_callback_LBD_Imp

Callback function used by SNOPT in a lower bounding problem (NLP) constructed by
relaxing the subproblem using the differentiable McCormick operators presented
by Khan2017. Inputs:
* `y::Vector{Float64}`: Point to evaluate in X
* `X::Vector{Interval{Float64}}`: Node corresponding to lower problem
* `opts`: Option type containing problem information
* `param`: Parameter object for generating implicit bounds
* `pmid`: Point at which parameter objects for implicit bounds are generated
Returns a tuple tuple (f_cv, c_cv, dfdx_cv, dcdx_cv, fail) where
* `f_cv::Float64`: Convex relaxation of the objective function
* `c_cv::Array{Float64,1}`: Convex relaxation of the constraint function
* `dfdx_cv::Array{Float64,1}`: Gradient of the convex relaxation of the
                               objective function
* `dcdx_cv::Array{Float64,2}`: Jacobian of the convex relaxation of the
                               constraint function
* `fail::Bool`: Flag used by SNOPT solver to indicate if SNOPT solver failed.
"""

function snopt_callback_LBD_Imp(y::Vector{Float64},
                                Y::Vector{Interval{Float64}},
                                opt,param,pmid)

    nx::Int64 = opt.Imp_nx
    np::Int64 = opt.numVar - nx
    l::Vector{Float64} = [Y[nx+i].lo for i=1:np]
    u::Vector{Float64} = [Y[nx+i].hi for i=1:np]
    pmid::Vector{Float64} = (l + u)/2.0

    # Evaluates the relaxation of objective and constraints
    # relaxation of function
    if opt.Imp_nCons>0
        f_val,c = impRelax_fg(opt.Imp_f,
                              opt.Imp_g,
                              opt.Imp_h,
                              opt.Imp_hj,
                              Y[1:nx],
                              Y[(nx+1):opt.numVar],
                              y,
                              pmid,
                              opt.solver.PSmcOpt,
                              param)
        dcdx::SparseMatrixCSC{Float64,Int64} = spzeros(length(opt.Imp_gL_loc)+length(opt.Imp_gU_loc),np)
    else
        f_val = impRelax_f(opt.Imp_f,
                           opt.Imp_h,
                           opt.Imp_hj,
                           Y[1:nx],
                           Y[(nx+1):opt.numVar],
                           y,
                           pmid,
                           opt.solver.PSmcOpt,
                           param)
        dcdx = spzeros(1,np)
    end

    # forms coefficients of linear relaxations
    if opt.Imp_nCons>0
        cx_ind1::Int64 = 1
        for i in opt.Imp_gL_loc
            for j=1:np
                if (c[i].cv_grad[j] != 0.0)
                    dcdx[cx_ind1,j] = c[i].cv_grad[j]
                end
            end
            cx_ind1 += 1
        end
        for i in opt.Imp_gU_loc
            for j=1:np
                if (c[i].cc_grad[j] != 0.0)
                    dcdx[cx_ind1,j] = -c[i].cc_grad[j]
                end
            end
            cx_ind1 += 1
        end
    end

    # forms rhs of linear relaxations
    if opt.Imp_nCons>0
        rhs::Vector{Float64} = zeros(Float64,length(opt.Imp_gL_loc)+length(opt.Imp_gU_loc))
    else
        rhs = zeros(Float64,1)
    end

    if opt.Imp_nCons>0
        cx_ind2::Int64 = 1
        for i in opt.Imp_gU_loc
            rhs[cx_ind2] = sum(pmid[:].*c[i].cv_grad[:])+opt.Imp_gU[i]-c[i].cv
            cx_ind2 += 1
        end
        for i in opt.Imp_gL_loc
            rhs[cx_ind2] = sum(-pmid[:].*c[i].cc_grad[:])-opt.Imp_gL[i]+c[i].cc
            cx_ind2 += 1
        end
    end

    fail::Bool = false
    return f_val.cv, rhs, f_val.cv_grad, dcdx, fail
end

function snopt_callback_LBD_Imp(y::Vector{Float64},
                                Y::Vector{MCInterval{Float64}},
                                opt,param,pmid)

    nx::Int64 = opt.solver.Imp_nx
    np::Int64 = opt.numVar - nx
    l::Vector{Float64} = [Y[nx+i].lo for i=1:np]
    u::Vector{Float64} = [Y[nx+i].hi for i=1:np]
    pmid::Vector{Float64} = (l + u)/2.0

    # Evaluates the relaxation of objective and constraints
    # relaxation of function

    if opt.Imp_nCons>0
        f_val,c = impRelax_fg(opt.Imp_f,
                              opt.Imp_g,
                              opt.Imp_h,
                              opt.Imp_hj,
                              Y[1:nx],
                              Y[(nx+1):opt.numVar],
                              y,
                              pmid,
                              opt.PSmcOpt,
                              param)
        dcdx::SparseMatrixCSC{Float64,Int64} = spzeros(length(opt.Imp_gL_loc)+length(opt.Imp_gU_loc),np)
    else
        f_val = impRelax_f(opt.Imp_f,
                           opt.Imp_h,
                           opt.Imp_hj,
                           Y[1:nx],
                           Y[(nx+1):opt.numVar],
                           y,
                           pmid,
                           opt.PSmcOpt,
                           param)
        dcdx = spzeros(1,np)
    end

    # forms coefficients of linear relaxations
    if opt.Imp_nCons>0
        cx_ind1::Int64 = 1
        for i in opt.Imp_gL_loc
            for j=1:np
                if (c[i].cv_grad[j] != 0.0)
                    dcdx[cx_ind1,j] = c[i].cv_grad[j]
                end
            end
            cx_ind1 += 1
        end
        for i in opt.Imp_gU_loc
            for j=1:np
                if (c[i].cc_grad[j] != 0.0)
                    dcdx[cx_ind1,j] = -c[i].cc_grad[j]
                end
            end
            cx_ind1 += 1
        end
    end

    # forms rhs of linear relaxations
    if opt.Imp_nCons>0
        rhs::Vector{Float64} = zeros(Float64,length(opt.Imp_gL_loc)+length(opt.Imp_gU_loc))
    else
        rhs = zeros(Float64,1)
    end

    if opt.Imp_nCons>0
        cx_ind2::Int64 = 1
        for i in opt.Imp_gU_loc
            rhs[cx_ind2] = sum(pmid[:].*c[i].cv_grad[:])+opt.Imp_gU[i]-c[i].cv
            cx_ind2 += 1
        end
        for i in opt.Imp_gL_loc
            rhs[cx_ind2] = sum(-pmid[:].*c[i].cc_grad[:])-opt.Imp_gL[i]+c[i].cc
            cx_ind2 += 1
        end
    end

    fail::Bool = false
    return f_val.cv, rhs, f_val.cv_grad, dcdx, fail
end

"""
    SNOPT_LBD_Imp

Solves a lower bounding problem (NLP) constructed by relaxing the subproblem
problem using the differentiable McCormick operators presented by Khan2017. The
relaxed problem is then solved using SNOPT and using the implicit bounding routine
of Stuber2014. SNOPT options and domain reduction settings are controlled in the
B&B algorithm using the global_options type. Inputs:
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
function SNOPT_LBD_Imp(Y::Vector{Interval{Float64}},
                       k::Int64,
                       pos::Int64,
                       opt,
                       UBD::Float64)

        options = Dict{String, Any}()
        options["Function precision"] = 1.0E-10
        options["Derivative option"] = 1
        options["Print file"] = 0
        options["Major print level"] = 0
        options["Minor print level"] = 0
        #options["Major optimality tolerance"] = 1e-6
        options["Verify level"] = 0
        nx::Int64 = opt[1].Imp_nx
        np::Int64 = opt[1].numVar - nx

        l::Vector{Float64} = [Y[nx+i].lo for i=1:np]
        u::Vector{Float64} = [Y[nx+i].hi for i=1:np]
        pmid::Vector{Float64} = (l + u)/2.0
        param = GenExpansionParams(opt[1].Imp_h,
                                   opt[1].Imp_hj,
                                   Y[1:nx],Y[(nx+1):opt[1].numVar],pmid,
                                   opt[1].solver.PSmcOpt)
        pnt::Vector{Float64}, val::Float64, status::Int64, mult::Vector{Float64} = Snopt.snopt(y::Vector{Float64} -> opt[2].fg_SNOPT_LBD_Imp(y,Y,param,pmid), pmid, l, u, options)
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
            for i=1:np
                if abs(mult[i])>opt[1].solver.dual_tol
                    if (pnt[i]-Y[i].lo)==(Y[i].hi-pnt[i])
                        mult_lo[i] = mult[i]
                        mult_hi[i] = mult[i]
                    elseif (pnt[i]-Y[i].lo)<(Y[i].hi-pnt[i])
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

function SNOPT_LBD_Imp(Y::Vector{MCInterval{Float64}},
                       k::Int64,
                       pos::Int64,
                       opt,
                       UBD::Float64)
        options = Dict{String, Any}()
        options["Function precision"] = 1.0E-10
        options["Derivative option"] = 1
        options["Print file"] = 0
        options["Major print level"] = 0
        options["Minor print level"] = 0
        #options["Major optimality tolerance"] = 1e-6
        options["Verify level"] = 0
        nx::Int64 = opt[1].solver.Implicit_Options.nx
        np::Int64 = opt[1].numVar - nx

        l::Vector{Float64} = [Y[nx+i].lo for i=1:np]
        u::Vector{Float64} = [Y[nx+i].hi for i=1:np]
        pmid::Vector{Float64} = (l + u)/2.0
        param = GenExpansionParams(opt[1].Imp_h,
                                   opt[1].Imp_hj,
                                   Y[1:nx],Y[(nx+1):opt[1].numVar],pmid,
                                   opt[1].solver.PSmcOpt)
        pnt::Vector{Float64}, val::Float64, status::Int64, mult::Vector{Float64} = Snopt.snopt(y::Vector{Float64} -> opt[2].fg_SNOPT_LBD_Imp(y,Y,param,pmid), pmid, l, u, options)
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
            for i=1:np
                if abs(mult[i])>opt[1].solver.dual_tol
                    if (pnt[i]-Y[i].lo)==(Y[i].hi-pnt[i])
                        mult_lo[i] = mult[i]
                        mult_hi[i] = mult[i]
                    elseif (pnt[i]-Y[i].lo)<(Y[i].hi-pnt[i])
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
