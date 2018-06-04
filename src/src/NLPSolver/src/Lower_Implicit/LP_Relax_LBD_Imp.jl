"""
    LP_Relax_LBD_Imp

Solves a lower bounding problem (LP) constructed by relaxing using nonsmooth
McCormick operators presented by McCormick1978, Mitsos2008. The relaxed problem
is then solved via an LP solvers. LP solver options and domain reduction
settings are controlled in the B&B algorithm using the global_options type. Inputs:
* `X::Vector{Interval{Float64}}`: Node over which to solve the lower problem
* `k::Int64`: Number of Iterations The B&B Algorithm has taken
* `pos::Int64`: Depth of Node in B&B Tree
* `opt`: Option type containing problem information
* `UBD::Float64`: Global upper bound for B&B algorithm

Returns a tuple `(val,pnt,feas,X,[])` where
* `val::Float64`: Lower bound calculated
* `pnt::Array{Float64,1}`: An array of length equal to X that gives the
                           optimal solution of the lower bound problem.
* `feas::Bool`: Returns true if the problem is feasible and false if it is
                infeasible
* `X::Vector{Interval{Float64}}`: Node over which the problem was solved after
                                  potential use of domain reduction techniques
* `[]`: The last element of the tuple is currently unused for this option.
"""
function LP_Relax_LBD_Imp(Y::Vector{Interval{Float64}},
                          k::Int64,
                          pos::Int64,
                          opt,
                          UBD::Float64)
        nx::Int64 = opt[1].Imp_nx
        np::Int64 = opt[1].Imp_np
        #println("ran lower implicit 1")
        try
            l::Vector{Float64} = [Y[nx+i].lo for i=1:np]
            u::Vector{Float64} = [Y[nx+i].hi for i=1:np]
            pmid::Vector{Float64} = (l + u)/2.0

            opt[1].solver.SubGradRefine && set_hybrid_box!(SVector{opt[1].numVar,Interval{Float64}}([Y[nx+i] for i=1:np]),
                                                           SVector{opt[1].numVar,Float64}(pmid),true)
            #println("pmid: $pmid")
            param = GenExpansionParams(opt[1].Imp_h, opt[1].Imp_hj,
                                         Y[1:nx],Y[(nx+1):(opt[1].numVar)],pmid,
                                         opt[1].solver.PSmcOpt)
            #=
            x_mc = param[end]
            p_mc::Vector{SMCg{np,Interval{Float64},Float64}} = [SMCg{np,Interval{Float64},Float64}(pmid[i],
                                                           pmid[i],
                                                           seed_g(Float64,i,np),
                                                           seed_g(Float64,i,np),
                                                           Y[nx+i],
                                                           false) for i=1:np]
            f::SMCg{np,Interval{Float64},Float64} = opt[1].Imp_f(x_mc[1:nx],p_mc)
            f_cv::Float64 = f.cv
            =#

            x_mc::Vector{HybridMC{np,Interval{Float64},Float64}} = param[end]
            p_mc::Vector{HybridMC{np,Interval{Float64},Float64}} = [HybridMC{np,Interval{Float64},Float64}(
                                                                    SMCg{np,Interval{Float64},Float64}(pmid[i],
                                                                                                       pmid[i],
                                                                                                       seed_g(Float64,i,np),
                                                                                                       seed_g(Float64,i,np),
                                                                                                       Y[nx+i],
                                                                                                       false)) for i=1:np]
            f::HybridMC{np,Interval{Float64},Float64} = opt[1].Imp_f(x_mc[1:nx],p_mc)
            f_cv::Float64 = cv(f)

            if opt[1].Imp_nCons>0
                #c::Vector{SMCg{np,Interval{Float64},Float64}} = opt[1].Imp_g(x_mc[1:nx],p_mc)
                c::Vector{HybridMC{np,Interval{Float64},Float64}} = opt[1].Imp_g(x_mc[1:nx],p_mc)
                dcdx::SparseMatrixCSC{Float64,Int64} = spzeros(length(opt[1].Imp_gL_Loc)+length(opt[1].Imp_gU_Loc),np)
            else
                dcdx = spzeros(1,np)
            end
            if opt[1].Imp_nCons>0
                cx_ind1 = 1
                for i in opt[1].Imp_gU_Loc
                    for j=1:np
                        if (c[i].cv_grad[j] != 0.0)
                            dcdx[cx_ind1,j] = c[i].cv_grad[j]
                        end
                    end
                    cx_ind1 += 1
                end
                for i in opt[1].Imp_gL_Loc
                    for j=1:np
                        if (c[i].cc_grad[j] != 0.0)
                            dcdx[cx_ind1,j] = -c[i].cc_grad[j]
                        end
                    end
                    cx_ind1 += 1
                end
            end
            if opt[1].Imp_nCons>0
                rhs::Vector{Float64} = zeros(Float64,length(opt[1].Imp_gL_Loc)+length(opt[1].Imp_gU_Loc))
            else
                rhs = zeros(Float64,1)
            end
            if opt[1].Imp_nCons>0
                cx_ind2 = 1
                for i in opt[1].Imp_gU_Loc
                    rhs[cx_ind2] = sum(pmid[:].*c[i].cv_grad[:])+opt[1].Imp_gU[i]-c[i].cv
                    cx_ind2 += 1
                end
                for i in opt[1].Imp_gL_Loc
                    rhs[cx_ind2] = sum(-pmid[:].*c[i].cc_grad[:])-opt[1].Imp_gL[i]+c[i].cc
                    cx_ind2 += 1
                end
            end
            model = buildlp([f.cv_grad[i] for i=1:np], dcdx, '<', rhs, l, u, opt[1].solver.LP_solver)
            result = solvelp(model)
            if (result.status == :Optimal)
                val::Float64 = max(f.Intv.lo,result.objval + f_cv - sum([pmid[i]*f.cv_grad[i] for i=1:np]))
                pnt::Vector{Float64} = vcat(mid.(Intv.(x_mc)),result.sol)
                feas::Bool = true
                mult::Vector{Float64} = result.attrs[:redcost]
                if (opt[1].Imp_nCons < 1)
                else
                    GInt::Vector{Interval{Float64}} = opt[1].Imp_g(Y[1:nx],Y[(nx+1):end])
                    cInt::Vector{Interval{Float64}} = vcat(GInt[opt[1].Imp_gU_Loc]-opt[1].Imp_gU[opt[1].Imp_gU_Loc],
                                                          -GInt[opt[1].Imp_gL_Loc]+opt[1].Imp_gL[opt[1].Imp_gL_Loc])
                    for i=1:length(cInt)
                        if (cInt[i].lo>0.0)
                            val = -Inf
                            pnt = pmid
                            feas = false
                            mult = pmid
                            break
                        end
                    end
                end
            elseif (result.status == :Infeasible)
                val = -Inf
                pnt = pmid
                feas = false
                mult = pmid
            else
                error("Solver error code $(result.status) in solver. Solution routine terminated.")
            end

            mult_lo::Vector{Float64} = [tol_eq(l[i],pnt[i],opt[1].solver.dual_tol) ? mult[i] : 0.0 for i=1:np]
            mult_hi::Vector{Float64} = [tol_eq(u[i],pnt[i],opt[1].solver.dual_tol) ? mult[i] : 0.0 for i=1:np]
            temp = Any[mult_lo,mult_hi,val]
            return val, pnt, feas, temp
        catch
            FInt::Interval = opt[1].Imp_f(Y[1:nx],Y[(nx+1):end])
            feas = true
            if (opt[1].Imp_nCons < 1)
            else
                GInt::Vector{Interval{Float64}} = opt[1].Imp_g(Y[1:nx],Y[(nx+1):end])
                cInt::Vector{Interval{Float64}} = vcat(GInt[opt[1].Imp_gU_Loc]-opt[1].Imp_gU[opt[1].Imp_gU_Loc],
                                                      -GInt[opt[1].Imp_gL_Loc]+opt[1].Imp_gL[opt[1].Imp_gL_Loc])
                for i=1:length(cInt)
                    if (cInt[i].lo>0.0)
                        feas = false
                        break
                    end
                end
            end
            mult_lo = [0.0 for i=1:np]
            mult_hi = [0.0 for i=1:np]
            val = FInt.lo
            pnt = mid.(Y)
            temp = Any[mult_lo,mult_hi,val]
            return val, pnt, feas, temp
        end
end

function LP_Relax_LBD_Imp(Y::Vector{MCInterval{Float64}},
                          k::Int64,
                          pos::Int64,
                          opt,
                          UBD::Float64)
        nx::Int64 = opt[1].Imp_nx
        np::Int64 = opt[1].numVar - nx
        try
            l::Vector{Float64} = [Y[nx+i].lo for i=1:np]
            u::Vector{Float64} = [Y[nx+i].hi for i=1:np]
            pmid::Vector{Float64} = (l + u)/2.0
            param = GenExpansionParams(opt[1].Imp_h, opt[1].Imp_hj,
                                         Y[1:nx],Y[(nx+1):(opt[1].numVar)],pmid,
                                         opt[1].solver.PSmcOpt)
            #=
            x_mc = param[end]
            p_mc::Vector{SMCg{np,MCInterval{Float64},Float64}} = [SMCg{np,MCInterval{Float64},Float64}(pmid[i],
                                                           pmid[i],
                                                           seed_g(Float64,i,np),
                                                           seed_g(Float64,i,np),
                                                           Y[nx+i],
                                                           false) for i=1:np]
            f::SMCg{np,MCInterval{Float64},Float64} = opt[1].Imp_f(x_mc[1:nx],p_mc)
            f_cv::Float64 = f.cv
            =#

            opt[1].solver.SubGradRefine && set_hybrid_box!([Y[nx+i] for i=1:np],pmid,true)
            x_mc::Vector{HybridMC{np,MCInterval{Float64},Float64}} = param[end]
            p_mc::Vector{HybridMC{np,MCInterval{Float64},Float64}} = [HybridMC{np,MCInterval{Float64},Float64}(
                                                                      SMCg{np,MCInterval{Float64},Float64}(pmid[i],
                                                                                                           pmid[i],
                                                                                                           seed_g(Float64,i,np),
                                                                                                           seed_g(Float64,i,np),
                                                                                                           Y[nx+i],
                                                                                                           false)) for i=1:np]
            f::HybridMC{np,MCInterval{Float64},Float64} = opt[1].Imp_f(x_mc[1:nx],p_mc)
            f_cv::Float64 = cv(f)

            if opt[1].Imp_nCons>0
                #c::Vector{SMCg{np,MCInterval{Float64},Float64}} = opt[1].Imp_g(x_mc[1:nx],p_mc)
                c::Vector{HybridMC{np,MCInterval{Float64},Float64}} = opt[1].Imp_g(x_mc[1:nx],p_mc)
                dcdx::SparseMatrixCSC{Float64,Int64} = spzeros(length(opt[1].Imp_gL_Loc)+length(opt[1].Imp_gU_Loc),np)
            else
                dcdx = spzeros(1,np)
            end
            if opt[1].Imp_nCons>0
                cx_ind1::Int64 = 1
                for i in opt[1].Imp_gL_Loc
                    for j=1:np
                        if (c[i].cv_grad[j] != 0.0)
                            dcdx[cx_ind1,j] = c[i].cv_grad[j]
                        end
                    end
                    cx_ind1 += 1
                end
                for i in opt[1].Imp_gU_Loc
                    for j=1:np
                        if (c[i].cc_grad[j] != 0.0)
                            dcdx[cx_ind1,j] = -c[i].cc_grad[j]
                        end
                    end
                    cx_ind1 += 1
                end
            end
            if opt[1].Imp_nCons>0
                rhs::Vector{Float64} = zeros(Float64,length(opt[1].Imp_gL_Loc)+length(opt[1].Imp_gU_Loc))
            else
                rhs = zeros(Float64,1)
            end
            if opt[1].Imp_nCons>0
                cx_ind2::Int64 = 1
                for i in opt[1].Imp_gU_Loc
                    rhs[cx_ind2] = sum(pmid[:].*c[i].cv_grad[:])+opt[1].Imp_gU[i]-c[i].cv
                    cx_ind2 += 1
                end
                for i in opt[1].Imp_gL_Loc
                    rhs[cx_ind2] = sum(-pmid[:].*c[i].cc_grad[:])-opt[1].Imp_gL[i]+c[i].cc
                    cx_ind2 += 1
                end
            end
            model = buildlp([f.cv_grad[i] for i=1:np], dcdx, '<', rhs, l, u, opt[1].solver.LP_solver)
            result = solvelp(model)
            if (result.status == :Optimal)
                val::Float64 = result.objval + f_cv - sum([pmid[i]*f.cv_grad[i] for i=1:np])
                pnt::Vector{Float64} = vcat(mid.(Intv.(x_mc)),result.sol)
                feas::Bool = true
                mult::Vector{Float64} = result.attrs[:redcost]
            elseif (result.status == :Infeasible)
                val = -Inf
                pnt = pmid
                feas = false
                mult = pmid
            else
                error("Solver error code $(result.status) in solver. Solution routine terminated.")
            end

            mult_lo::Vector{Float64} = [tol_eq(l[i],pnt[i],opt[1].solver.dual_tol) ? mult[i] : 0.0 for i=1:np]
            mult_hi::Vector{Float64} = [tol_eq(u[i],pnt[i],opt[1].solver.dual_tol) ? mult[i] : 0.0 for i=1:np]
            temp = Any[mult_lo,mult_hi,val]
            return val, pnt, feas, temp
        catch
            FInt::MCInterval{Float64} = opt[1].Imp_f(Y[1:nx],Y[(nx+1):end])
            feas = true
            if (opt[1].Imp_nCons < 1)
            else
                GInt::Vector{MCInterval{Float64}} = opt[1].Imp_g(Y[1:nx],Y[(nx+1):end])
                cInt::Vector{MCInterval{Float64}} = vcat(GInt[opt[1].Imp_gU_Loc]-opt[1].Imp_gU[opt[1].Imp_gU_Loc],
                                                        -GInt[opt[1].Imp_gL_Loc]+opt[1].Imp_gL[opt[1].Imp_gL_Loc])
                for i=1:length(cInt)
                    if (cInt[i].lo>0.0)
                        feas = false
                        break
                    end
                end
            end
            mult_lo = [0.0 for i=1:np]
            mult_hi = [0.0 for i=1:np]
            val = FInt.lo
            pnt = mid.(Y)
            temp = Any[mult_lo,mult_hi,val]
            return val, pnt, feas, temp
        end
end
