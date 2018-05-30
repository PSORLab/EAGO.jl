"""
    LP_Relax_LBD

Solves a lower bounding problem (LP) constructed by relaxing using nonsmooth
McCormick operators presented by McCormick1978, Mitsos2008. The relaxed problem
is then solved via an LP solvers. LP solver options and domain reduction
settings are controlled in the B&B algorithm using the global_options type. Inputs:
* `X::Vector{Interval{Float64}}`: The interval bounds of the relaxated LBD problem.
* `k::Int64`: Iteration number of the B&B problem.
* `pos::Int64`: Position of the node in the B&B tree.
* `opt::Any`: Option storage for problem
* `UBD::Float64`: Upper bound
This returns a tuple `(val,pnt,feas,X,[mult_lo,mult_hi,val])` where
* `val::Float64`: Lower bound calculated
* `pnt::Vector{Float64}`: An array of length equal to X that gives the optimal
                          solution of the lower bound problem.
* `feas::Bool`: Returns true if the problem is feasible and false if it is infeasible
* `X::Vector{Interval}`: Node over which the problem was solved after
                          potential use of domain reduction techniques
* `[mult_lo,mult_hi,val]`: Info for OBBT.
"""
function LP_Relax_LBD(X::Vector{Interval{Float64}},
                      k::Int64,
                      pos::Int64,
                      opt::Any,
                      UBD::Float64)

        l::Vector{Float64} = [X[i].lo for i=1:opt[1].numVar]
        u::Vector{Float64} = [X[i].hi for i=1:opt[1].numVar]
        x0::Vector{Float64} = (l + u)/2.0
        x_mc::Vector{SMCg{opt[1].numVar,Interval{Float64},Float64}} = [SMCg{opt[1].numVar,Interval{Float64},Float64}(x0[i],
                                                                            x0[i],
                                                                            seed_g(Float64,i,opt[1].numVar),
                                                                            seed_g(Float64,i,opt[1].numVar),
                                                                            X[i],
                                                                            false) for i=1:opt[1].numVar]

        # relaxation of function
        f::SMCg{opt[1].numVar,Interval{Float64},Float64} = opt[1].f(x_mc)
        f_cv::Float64 = f.cv

        if opt[1].numConstr>0
            c::Vector{SMCg{opt[1].numVar,Interval{Float64},Float64}} = opt[1].g(x_mc)
            #println("c: $c")
            dcdx::SparseMatrixCSC{Float64,Int64} = spzeros(length(opt[1].gL_loc)+length(opt[1].gU_loc),opt[1].numVar)
        else
            dcdx = spzeros(1,opt[1].numVar)
        end

        # forms coefficients of linear relaxations
        if opt[1].numConstr>0
            cx_ind1::Int64 = 1
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
            rhs::Vector{Float64} = zeros(Float64,length(opt[1].gL_loc)+length(opt[1].gU_loc))
        else
            rhs = zeros(Float64,1)
        end

        if opt[1].numConstr>0
            cx_ind2::Int64 = 1
            for i in opt[1].gU_loc
                rhs[cx_ind2] = sum(x0[:].*c[i].cv_grad[:])+min(c[i].Intv.hi,opt[1].gU[i])-c[i].cv
                cx_ind2 += 1
            end
            for i in opt[1].gL_loc
                rhs[cx_ind2] = sum(-x0[:].*c[i].cc_grad[:])-max(c[i].Intv.lo,opt[1].gL[i])+c[i].cc
                cx_ind2 += 1
            end
        end

        model = buildlp([f.cv_grad[i] for i=1:opt[1].numVar], dcdx, '<', rhs, l, u, opt[1].solver.LP_solver)
        result = solvelp(model)

        # Unpacks the results from the LP solver
        if (result.status == :Optimal)
            val::Float64 = max(f.Intv.lo,result.objval + f_cv - sum([x0[i]*f.cv_grad[i] for i=1:opt[1].numVar]))
            pnt::Vector{Float64} = result.sol
            feas = true
            mult::Vector{Float64} = result.attrs[:redcost]

            if (opt[1].numConstr > 0)
              GInt = opt[1].g(X)
              cInt = vcat(GInt[opt[1].gU_loc]-opt[1].gU[opt[1].gU_loc],
                                    -GInt[opt[1].gL_loc]+opt[1].gL[opt[1].gL_loc])
              for i=1:opt[1].gexp
                if (cInt[i].lo>0.0)
                    val = -Inf
                    pnt = x0
                    feas = false
                    mult = x0
                  break
                end
              end
            end

        elseif (result.status == :Infeasible)
            val = -Inf
            pnt = x0
            feas = false
            mult = x0
        else
            error("Solver error code $(result.status) in solver. Solution routine terminated.")
        end

        # Formats Duality-Based Multipliers
        mult_lo = [tol_eq(X[i].lo,pnt[i],opt[1].solver.dual_tol) ? mult[i] : 0.0 for i=1:opt[1].numVar]
        mult_hi = [tol_eq(X[i].hi,pnt[i],opt[1].solver.dual_tol) ? mult[i] : 0.0 for i=1:opt[1].numVar]
        temp = Any[mult_lo,mult_hi,val]

    #println("END LOWER BOUNDING PROBLEM!")
    return val, pnt, feas, temp
end
function LP_Relax_LBD(X::Vector{MCInterval{Float64}},
                      k::Int64,
                      pos::Int64,
                      opt::Any,
                      UBD::Float64)

        l::Vector{Float64} = [X[i].lo for i=1:opt[1].numVar]
        u::Vector{Float64} = [X[i].hi for i=1:opt[1].numVar]
        x0::Vector{Float64} = (l + u)/2.0
        x_mc::Vector{SMCg{opt[1].numVar,MCInterval{Float64},Float64}} = [SMCg{opt[1].numVar,MCInterval{Float64},Float64}(x0[i],
                                                                            x0[i],
                                                                            seed_g(Float64,i,opt[1].numVar),
                                                                            seed_g(Float64,i,opt[1].numVar),
                                                                            X[i],
                                                                            false) for i=1:opt[1].numVar]
        # relaxation of function
        f::SMCg{opt[1].numVar,MCInterval{Float64},Float64} = opt[1].f(x_mc)
        f_cv::Float64 = f.cv

        if opt[1].numConstr>0
            c::Vector{SMCg{opt[1].numVar,MCInterval{Float64},Float64}} = opt[1].g(x_mc)
            dcdx::SparseMatrixCSC{Float64,Int64} = spzeros(length(opt[1].gL_loc)+length(opt[1].gU_loc),opt[1].numVar)
        else
            dcdx = spzeros(1,opt[1].numVar)
        end

        # forms coefficients of linear relaxations
        if opt[1].numConstr>0
            cx_ind1::Int64 = 1
            for i in opt[1].gU_loc
                for j=1:opt[1].numVar
                    if (c[i].cv_grad[j] != 0.0)
                        dcdx[cx_ind1,j] = c[i].cv_grad[j]
                    end
                end
                cx_ind1 += 1
            end
            for i in opt[1].gL_loc
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
            rhs::Vector{Float64} = zeros(Float64,length(opt[1].gL_loc)+length(opt[1].gU_loc))
        else
            rhs = zeros(Float64,1)
        end

        if opt[1].numConstr>0
            cx_ind2::Int64 = 1
            for i in opt[1].gU_loc
                rhs[cx_ind2] = sum(x0[:].*c[i].cv_grad[:])+min(c[i].Intv.hi,opt[1].gU[i])-c[i].cv
                cx_ind2 += 1
            end
            for i in opt[1].gL_loc
                rhs[cx_ind2] = sum(-x0[:].*c[i].cc_grad[:])-max(c[i].Intv.lo,opt[1].gL[i])+c[i].cc
                cx_ind2 += 1
            end
        end

        model = buildlp([f.cv_grad[i] for i=1:opt[1].numVar], dcdx, '<', rhs, l, u, opt[1].solver.LP_solver)
        result = solvelp(model)
        # Unpacks the results from the LP solver
        if (result.status == :Optimal)
            val::Float64 = max(f.Intv.lo,result.objval + f_cv - sum([x0[i]*f.cv_grad[i] for i=1:opt[1].numVar]))
            pnt::Vector{Float64} = result.sol
            feas = true
            mult::Vector{Float64} = result.attrs[:redcost]


                        if (opt[1].numConstr > 0)
                          GInt = opt[1].g(X)
                          cInt = vcat(GInt[opt[1].gU_loc]-opt[1].gU[opt[1].gU_loc],
                                                -GInt[opt[1].gL_loc]+opt[1].gL[opt[1].gL_loc])
                          for i=1:opt[1].gexp
                            if (cInt[i].lo>0.0)
                                val = -Inf
                                pnt = x0
                                feas = false
                                mult = x0
                              break
                            end
                          end
                        end

        elseif (result.status == :Infeasible)
            val = -Inf
            pnt = x0
            feas = false
            mult = x0
        else
            error("Solver error code $(result.status) in solver. Solution routine terminated.")
        end

        # Formats Duality-Based Multipliers
        mult_lo = [tol_eq(X[i].lo,pnt[i],opt[1].solver.dual_tol) ? mult[i] : 0.0 for i=1:opt[1].numVar]
        mult_hi = [tol_eq(X[i].hi,pnt[i],opt[1].solver.dual_tol) ? mult[i] : 0.0 for i=1:opt[1].numVar]
        temp = Any[mult_lo,mult_hi,val]

    #println("END LOWER BOUNDING PROBLEM!")
    return val, pnt, feas, temp
end
