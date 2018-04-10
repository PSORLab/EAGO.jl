"""
    STD_LP_Probe!(X::Vector{Interval{Float64}},opt,UBD::Float64)

Performs standard probing using linear underestimator forms via
McCormick relaxations, sparse calculations, and contracts `X::Vector{Interval{Float64}}` in place.
The `opt` object has the fields `.numVar`, `.numConstr`, `.f`, `.g`, `.gL`, `.gU`,
`.gL_Loc`, `.gU_Loc` and `solver.LP_solver`. `f(x)` is the objective, `g(x)` is
the constraint function, `numVar` is the number of decision variables, `numConstr`
is the number of constraints, `gL` is the lower bound, `gU` is the upper bound,
`gL_Loc` is an index to indicate the lower bound is finite, and `gU_Loc` is an index
that indicates the upper bound is finite. The upper bound is `UBD::Float64`.
"""
function STD_LP_Probe!(X::Vector{Interval{Float64}},opt,UBD::Float64)

    mu_temp::Int64 = copy(EAGOSmoothMcCormickGrad.MC_param.mu)
    set_diff_relax(0)

    # constructs LP relaxation
    Xlo::Vector{Float64} = [X[i].lo for i=1:opt[1].numVar]
    Xhi::Vector{Float64} = [X[i].hi for i=1:opt[1].numVar]
    x0::Vector{Float64} = (Xlo + Xhi)/2.0
    x_SMC::Vector{SMCg{opt[1].numVar,Float64}} = [SMCg{opt[1].numVar,Float64}(x0[i],
                                                                              x0[i],
                                                                              seed_g(Float64,i,opt[1].numVar),
                                                                              seed_g(Float64,i,opt[1].numVar),
                                                                              X[i],
                                                                              false,
                                                                              X,
                                                                              x0) for i=1:opt[1].numVar]
    # probes upper bound
    f_mc::SMCg{opt[1].numVar,Float64} = opt[1].f(x_SMC)
    f_cv::Float64 = f_mc.cv
    if opt[1].numConstr>0
        c::Vector{SMCg{opt[1].numVar,Float64}} = opt[1].g(x_SMC)
        dcdx::SparseMatrixCSC{Float64,Int64} = spzeros(length(opt[1].gL_loc)+length(opt[1].gU_loc)+1,opt[1].numVar)
    else
        dcdx = spzeros(1,opt[1].numVar)
    end

    dcdx[1,:] = f_mc.cv_grad
    if opt[1].numConstr>0
        cx_ind1::Int64 = 2
        for i in opt[1].gL_loc
            for j=1:opt[1].numVar
                if (c[i].cv_grad[j] != 0.0)
                    dcdx[cx_ind1,j] = c[i].cv_grad[j]
                end
            end
            cx_ind1 += 1
        end
        for i in opt[1].gU_loc
            for j=1:opt[1].numVar
                if (c[i].cc_grad[j] != 0.0)
                    dcdx[cx_ind1,j] = -c[i].cc_grad[j]
                end
            end
            cx_ind1 += 1
        end
    end

    # forms rhs of linear relaxations
    if opt[1].numConstr>0
        rhs::Vector{Float64} = zeros(Float64,length(opt[1].gL_loc)+length(opt[1].gU_loc)+1)
    else
        rhs = zeros(Float64,1)
    end
    rhs[1] = UBD + f_mc.cv - sum(x0.*f_mc.cv_grad)

    if opt[1].numConstr>0
        cx_ind2::Int64 = 2
        for i in opt[1].gU_loc
            rhs[cx_ind2] = sum(x0[:].*c[i].cv_grad[:])+opt[1].gU[i]-c[i].cv
            cx_ind2 += 1
        end
        for i in opt[1].gL_loc
            rhs[cx_ind2] = sum(-x0[:].*c[i].cc_grad[:])-opt[1].gL[i]+c[i].cc
            cx_ind2 += 1
        end
    end

    if (opt[1].numConstr>0)
        temp_model = buildlp([f_mc.cv_grad[i] for i=1:opt[1].numVar], dcdx, '<', rhs, Xlo, Xhi, opt[1].solver.LP_solver)
    else
        temp_model = buildlp([f_mc.cv_grad[i] for i=1:opt[1].numVar], zeros(opt[1].numVar,opt[1].numVar), '<', zeros(opt[1].numVar), Xlo, Xhi, opt[1].solver.LP_solver)
    end

    for i=1:opt[1].numVar

        # probes upper bound
        MathProgBase.setvarLB!(temp_model,[i == j ? Xhi[i] : Xlo[i] for j=1:opt[1].numVar])
        MathProgBase.setvarUB!(temp_model,Xhi)
        result = solvelp(temp_model)
        if (result.status == :Optimal)
            val::Float64 = result.objval + f_cv - sum([x0[i]*f_mc.cv_grad[i] for i=1:opt[1].numVar])
            pnt::Vector{Float64} = result.sol
            mult::Vector{Float64} = result.attrs[:redcost]
            mult_lo::Vector{Float64} = [tol_eq(X[i].lo,pnt[i],opt[1].solver.dual_tol) ? mult[i] : 0.0 for i=1:opt[1].numVar]
            mult_hi::Vector{Float64} = [tol_eq(X[i].hi,pnt[i],opt[1].solver.dual_tol) ? mult[i] : 0.0 for i=1:opt[1].numVar]
            Variable_DR!(X,mult_lo,mult_hi,val,UBD)
        end
        # probes lower bound
        MathProgBase.setvarLB!(temp_model,Xlo)
        MathProgBase.setvarUB!(temp_model,[i == j ? Xhi[i] : Xlo[i] for j=1:opt[1].numVar])
        result = solvelp(temp_model)
        if (result.status == :Optimal)
            val = result.objval + f_cv - sum([x0[i]*f_mc.cv_grad[i] for i=1:opt[1].numVar])
            pnt = result.sol
            mult = result.attrs[:redcost]
            mult_lo = [tol_eq(X[i].lo,pnt[i],opt[1].solver.dual_tol) ? mult[i] : 0.0 for i=1:opt[1].numVar]
            mult_hi = [tol_eq(X[i].hi,pnt[i],opt[1].solver.dual_tol) ? mult[i] : 0.0 for i=1:opt[1].numVar]
            Variable_DR!(X,mult_lo,mult_hi,val,UBD)
        end
    end

    # reset McCormick relaxation to original value
    set_diff_relax(mu_temp)

    # stores outputs
    X[:] = X
    return true
end

"""
    Imp_LP_Probe!(X::Vector{Interval{Float64}},opt,UBD::Float64)

Performs standard probing using linear underestimator forms via
McCormick relaxations, sparse calculations, and contracts `X::Vector{Interval{Float64}}`
in place only tightening the following intervals in the vector `X[(nx+1):(nx+np)]`.
The `opt` object has the fields `.numVar`, `.numConstr`, `.f`, `.g`, `.gL`, `.gU`,
`.gL_Loc`, `.gU_Loc` and `solver.LP_solver`. `f(x)` is the objective, `g(x)` is
the constraint function, `numVar` is the number of decision variables, `numConstr`
is the number of constraints, `gL` is the lower bound, `gU` is the upper bound,
`gL_Loc` is an index to indicate the lower bound is finite, and `gU_Loc` is an index
that indicates the upper bound is finite. The upper bound is `UBD::Float64`.
"""
function Imp_LP_Probe!(X::Vector{Interval{Float64}},opt,UBD::Float64)

    # constructs LP relaxation
    Xlo::Vector{Float64} = [X[i].lo for i=1:opt[1].numVar]
    Xhi::Vector{Float64} = [X[i].hi for i=1:opt[1].numVar]
    x0::Vector{Float64} = (Xlo + Xhi)/2.0
    x_SMC::Vector{SMCg{opt[1].numVar,Float64}} = [SMCg{opt[1].numVar,Float64}(x0[i],
                                                                              x0[i],
                                                                              seed_g(Float64,i,opt[1].numVar),
                                                                              seed_g(Float64,i,opt[1].numVar),
                                                                              X[i],
                                                                              false,
                                                                              X,
                                                                              x0) for i=1:opt[1].numVar]
    # probes upper bound
    f_mc::SMCg{opt[1].numVar,Float64} = opt[1].f(x_SMC)
    f_cv::Float64 = f_mc.cv
    c::Vector{SMCg{opt[1].numVar,Float64}} = opt[1].g(x_SMC)
    # forms coefficients of linear relaxations
    dcdx::SparseMatrixCSC{Float64,Int64} = spzeros(length(opt[1].gL_loc)+length(opt[1].gU_loc)+1,opt[1].numVar)
    dcdx[1,:] = f_mc.cv_grad
    cx_ind1::Int64 = 2
    for i in opt[1].gL_loc
        for j=1:opt[1].numVar
            if (c[i].cv_grad[j] != 0.0)
                dcdx[cx_ind1,j] = c[i].cv_grad[j]
            end
        end
        cx_ind1 += 1
    end
    for i in opt[1].gU_loc
        for j=1:opt[1].numVar
            if (c[i].cc_grad[j] != 0.0)
                dcdx[cx_ind1,j] = -c[i].cc_grad[j]
            end
        end
        cx_ind1 += 1
    end

    # forms rhs of linear relaxations
    rhs::Vector{Float64} = zeros(Float64,length(opt[1].gL_loc)+length(opt[1].gU_loc)+1)
    rhs[1] = UBD + f_mc.cv - sum(x0.*f_mc.cv_grad)
    cx_ind2::Int64 = 2
    for i in opt[1].gU_loc
        rhs[cx_ind2] = sum(x0[:].*c[i].cv_grad[:])+opt[1].gU[i]-c[i].cv
        cx_ind2 += 1
    end
    for i in opt[1].gL_loc
        rhs[cx_ind2] = sum(-x0[:].*c[i].cc_grad[:])-opt[1].gL[i]+c[i].cc
        cx_ind2 += 1
    end

    if (opt[1].numConstr>0)
        temp_model = buildlp([f_mc.cv_grad[i] for i=1:opt[1].numVar], dcdx, '<', rhs, Xlo, Xhi, opt[1].solver.LP_solver)
    else
        temp_model = buildlp([f_mc.cv_grad[i] for i=1:opt[1].numVar], zeros(opt[1].numVar,opt[1].numVar), '<', zeros(opt[1].numVar), Xlo, Xhi, opt[1].solver.LP_solver)
    end

    for i=1:opt[1].numVar

        # probes upper bound
        MathProgBase.setvarLB!(temp_model,[i == j ? Xhi[i] : Xlo[i] for j=1:opt[1].numVar])
        MathProgBase.setvarUB!(temp_model,Xhi)
        result = solvelp(temp_model)
        if (result.status == :Optimal)
            val::Float64 = result.objval + f_cv - sum([x0[i]*f_mc.cv_grad[i] for i=1:opt[1].numVar])
            pnt::Vector{Float64} = result.sol
            mult::Vector{Float64} = result.attrs[:redcost]
            mult_lo::Vector{Float64} = [tol_eq(X[i].lo,pnt[i],1E-12) ? mult[i] : 0.0 for i=1:opt[1].numVar]
            mult_hi::Vector{Float64} = [tol_eq(X[i].hi,pnt[i],1E-12) ? mult[i] : 0.0 for i=1:opt[1].numVar]
            Variable_DR!(X,mult_lo,mult_hi,val,UBD)
        end
        # probes lower bound
        MathProgBase.setvarLB!(temp_model,Xlo)
        MathProgBase.setvarUB!(temp_model,[i == j ? Xhi[i] : Xlo[i] for j=1:opt[1].numVar])
        result = solvelp(temp_model)
        if (result.status == :Optimal)
            val = result.objval + f_cv - sum([x0[i]*f_mc.cv_grad[i] for i=1:opt[1].numVar])
            pnt = result.sol
            mult = result.attrs[:redcost]
            mult_lo = [tol_eq(X[i].lo,pnt[i],1E-12) ? mult[i] : 0.0 for i=1:opt[1].numVar]
            mult_hi = [tol_eq(X[i].hi,pnt[i],1E-12) ? mult[i] : 0.0 for i=1:opt[1].numVar]
            Variable_DR!(X,mult_lo,mult_hi,val,UBD)
        end
    end
end
