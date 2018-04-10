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
* `UBD::Float64: Global upper bound for B&B algorithm

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
        #println("     -----     Lower LP Imp Start     -----     ")
        nx::Int64 = opt[1].solver.Implicit_Options.nx
        np::Int64 = opt[1].numVar - nx
        try
            #println("ran try!")
            l::Vector{Float64} = [Y[nx+i].lo for i=1:np]
            u::Vector{Float64} = [Y[nx+i].hi for i=1:np]
            pmid::Vector{Float64} = (l + u)/2.0

            #println("ran to me 2:")
            if (opt[1].solver.Implicit_Options.Inplace)
                #println(" LP param gen in place")
                param = IndGenExpansionParams(opt[1].solver.Implicit_Options.h,
                                         opt[1].solver.Implicit_Options.hj,
                                         Y[1:nx],Y[(nx+1):opt[1].numVar],pmid,
                                         opt[1].solver.Implicit_Options.opts)
                #println(" LP param gen in place finish")
            else
                #println(" LP param gen out place")
                param = GenExpansionParams(opt[1].solver.Implicit_Options.h,
                                       opt[1].solver.Implicit_Options.hj,
                                       Y[1:nx],Y[(nx+1):opt[1].numVar],pmid,
                                       opt[1].solver.Implicit_Options.opts)
                #println(" LP param gen out place finish")
            end
            x_mc = param[end]
            #println("x_mc: $x_mc")
            p_mc::Vector{SMCg{np,Interval{Float64},Float64}} = [SMCg{np,Interval{Float64},Float64}(pmid[i],
                                                           pmid[i],
                                                           seed_g(Float64,i,np),
                                                           seed_g(Float64,i,np),
                                                           Y[nx+i],
                                                           false,
                                                           Y[(nx+1):(nx+np)],
                                                           pmid) for i=1:np]
            #println("p_mc: $p_mc")
            # relaxation of function
            #println("ran to me 3:")
            #println("opt[1].solver.Implicit_Options.f(x_mc[1:nx],p_mc): $(opt[1].solver.Implicit_Options.f(x_mc[1:nx],p_mc))")
            f::SMCg{np,Interval{Float64},Float64} = opt[1].solver.Implicit_Options.f(x_mc[1:nx],p_mc)
            f_cv::Float64 = f.cv
            #println("ran to me 3a:")
            if opt[1].solver.Implicit_Options.numConstr>0
                c::Vector{SMCg{np,Interval{Float64},Float64}} = opt[1].solver.Implicit_Options.g(x_mc[1:nx],p_mc)
                dcdx::SparseMatrixCSC{Float64,Int64} = spzeros(length(opt[1].solver.Implicit_Options.gL_loc)+length(opt[1].solver.Implicit_Options.gU_loc),np)
            else
                dcdx = spzeros(1,np)
            end
            #println("ran to me 4:")
            # forms coefficients of linear relaxations
            if opt[1].solver.Implicit_Options.numConstr>0
                cx_ind1::Int64 = 1
                for i in opt[1].solver.Implicit_Options.gL_loc
                    for j=1:np
                        if (c[i].cv_grad[j] != 0.0)
                            dcdx[cx_ind1,j] = c[i].cv_grad[j]
                        end
                    end
                    cx_ind1 += 1
                end
                for i in opt[1].solver.Implicit_Options.gU_loc
                    for j=1:np
                        if (c[i].cc_grad[j] != 0.0)
                            dcdx[cx_ind1,j] = -c[i].cc_grad[j]
                        end
                    end
                    cx_ind1 += 1
                end
            end
            #println("ran to me 5:")
            # forms rhs of linear relaxations
            if opt[1].solver.Implicit_Options.numConstr>0
                rhs::Vector{Float64} = zeros(Float64,length(opt[1].solver.Implicit_Options.gL_loc)+length(opt[1].solver.Implicit_Options.gU_loc))
            else
                rhs = zeros(Float64,1)
            end
            #println("ran to me 6:")
            if opt[1].solver.Implicit_Options.numConstr>0
                cx_ind2::Int64 = 1
                for i in opt[1].solver.Implicit_Options.gU_loc
                    rhs[cx_ind2] = sum(pmid[:].*c[i].cv_grad[:])+opt[1].solver.Implicit_Options.gU[i]-c[i].cv
                    cx_ind2 += 1
                end
                for i in opt[1].solver.Implicit_Options.gL_loc
                    rhs[cx_ind2] = sum(-pmid[:].*c[i].cc_grad[:])-opt[1].solver.Implicit_Options.gL[i]+c[i].cc
                    cx_ind2 += 1
                end
            end
            #println("ran to me 7:")
            #println("f.cv_grad: $(f.cv_grad)")
            #println("l: $(l)")
            #println("u: $(u)")
            #println("dcdx: $(dcdx)")
            model = buildlp([f.cv_grad[i] for i=1:np], dcdx, '<', rhs, l, u, opt[1].solver.LP_solver)
            result = solvelp(model)
            #println("result: $result.status")
            # Unpacks the results from the LP solver
            if (result.status == :Optimal)
                #println("ran me 8b")
                val::Float64 = result.objval + f_cv - sum([pmid[i]*f.cv_grad[i] for i=1:np])
                pnt::Vector{Float64} = vcat(mid.(Intv.(x_mc)),result.sol)
                feas = true
                mult::Vector{Float64} = result.attrs[:redcost]
                #println("val: $val")
                #println("pnt: $pnt")
                #println("feas: $feas")
                #println("mult: $mult")
                #println("ran me 8f")
            elseif (result.status == :Infeasible)
                #println("ran me 8b1")
                val = -Inf
                pnt = pmid
                feas = false
                mult = pmid
                #println("ran me 8f1")
            else
                #println("ran me 8b2")
                error("Solver error code $(result.status) in solver. Solution routine terminated.")
                #println("ran me 8b2")
            end

            mult_lo = [tol_eq(l[i],pnt[i],opt[1].solver.dual_tol) ? mult[i] : 0.0 for i=1:np]
            mult_hi = [tol_eq(u[i],pnt[i],opt[1].solver.dual_tol) ? mult[i] : 0.0 for i=1:np]
            #println("mult_lo: $mult_lo")
            #println("mult_hi: $mult_hi")
            temp = Any[mult_lo,mult_hi,val]
            #println("     -----     Lower LP Try Imp End     -----     ")
            return val, pnt, feas, temp
        catch
            # solve optimization problem via interval extension
            #println("ran catch!")
            FInt::Interval = opt[1].solver.Implicit_Options.f(Y[1:nx],Y[(nx+1):end])
            feas = true
            if (opt[1].solver.Implicit_Options.numConstr < 1)
            else
                GInt::Vector{Interval{Float64}} = opt[1].solver.Implicit_Options.g(Y[1:nx],Y[(nx+1):end])
                cInt::Vector{Interval{Float64}} = vcat(GInt[opt[1].solver.Implicit_Options.gU_loc]-opt[1].solver.Implicit_Options.gU[opt[1].solver.Implicit_Options.gU_loc],
                                                      -GInt[opt[1].solver.Implicit_Options.gL_loc]+opt[1].solver.Implicit_Options.gL[opt[1].solver.Implicit_Options.gL_loc])
                # Update To sparse later
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
            #println("     -----     Lower LP Catch Imp End     -----     ")
            return val, pnt, feas, temp
        end
end

function LP_Relax_LBD_Imp(Y::Vector{MCInterval{Float64}},
                          k::Int64,
                          pos::Int64,
                          opt,
                          UBD::Float64)
        #println("     -----     Lower LP Imp Start     -----     ")
        nx::Int64 = opt[1].solver.Implicit_Options.nx
        np::Int64 = opt[1].numVar - nx
        try
            #println("ran try!")
            l::Vector{Float64} = [Y[nx+i].lo for i=1:np]
            u::Vector{Float64} = [Y[nx+i].hi for i=1:np]
            pmid::Vector{Float64} = (l + u)/2.0

            #println("ran to me 2:")
            if (opt[1].solver.Implicit_Options.Inplace)
                #println(" LP param gen in place")
                param = IndGenExpansionParams(opt[1].solver.Implicit_Options.h,
                                         opt[1].solver.Implicit_Options.hj,
                                         Y[1:nx],Y[(nx+1):opt[1].numVar],pmid,
                                         opt[1].solver.Implicit_Options.opts)
                #println(" LP param gen in place finish")
            else
                #println(" LP param gen out place")
                param = GenExpansionParams(opt[1].solver.Implicit_Options.h,
                                       opt[1].solver.Implicit_Options.hj,
                                       Y[1:nx],Y[(nx+1):opt[1].numVar],pmid,
                                       opt[1].solver.Implicit_Options.opts)
                #println(" LP param gen out place finish")
            end
            x_mc = param[end]
            #println("x_mc: $x_mc")
            p_mc::Vector{SMCg{np,MCInterval{Float64},Float64}} = [SMCg{np,MCInterval{Float64},Float64}(pmid[i],
                                                           pmid[i],
                                                           seed_g(Float64,i,np),
                                                           seed_g(Float64,i,np),
                                                           Y[nx+i],
                                                           false,
                                                           Y[(nx+1):(nx+np)],
                                                           pmid) for i=1:np]
            #println("p_mc: $p_mc")
            # relaxation of function
            #println("ran to me 3:")
            #println("opt[1].solver.Implicit_Options.f(x_mc[1:nx],p_mc): $(opt[1].solver.Implicit_Options.f(x_mc[1:nx],p_mc))")
            f::SMCg{np,MCInterval{Float64},Float64} = opt[1].solver.Implicit_Options.f(x_mc[1:nx],p_mc)
            f_cv::Float64 = f.cv
            #println("ran to me 3a:")
            if opt[1].solver.Implicit_Options.numConstr>0
                c::Vector{SMCg{np,MCInterval{Float64},Float64}} = opt[1].solver.Implicit_Options.g(x_mc[1:nx],p_mc)
                dcdx::SparseMatrixCSC{Float64,Int64} = spzeros(length(opt[1].solver.Implicit_Options.gL_loc)+length(opt[1].solver.Implicit_Options.gU_loc),np)
            else
                dcdx = spzeros(1,np)
            end
            #println("ran to me 4:")
            # forms coefficients of linear relaxations
            if opt[1].solver.Implicit_Options.numConstr>0
                cx_ind1::Int64 = 1
                for i in opt[1].solver.Implicit_Options.gL_loc
                    for j=1:np
                        if (c[i].cv_grad[j] != 0.0)
                            dcdx[cx_ind1,j] = c[i].cv_grad[j]
                        end
                    end
                    cx_ind1 += 1
                end
                for i in opt[1].solver.Implicit_Options.gU_loc
                    for j=1:np
                        if (c[i].cc_grad[j] != 0.0)
                            dcdx[cx_ind1,j] = -c[i].cc_grad[j]
                        end
                    end
                    cx_ind1 += 1
                end
            end
            #println("ran to me 5:")
            # forms rhs of linear relaxations
            if opt[1].solver.Implicit_Options.numConstr>0
                rhs::Vector{Float64} = zeros(Float64,length(opt[1].solver.Implicit_Options.gL_loc)+length(opt[1].solver.Implicit_Options.gU_loc))
            else
                rhs = zeros(Float64,1)
            end
            #println("ran to me 6:")
            if opt[1].solver.Implicit_Options.numConstr>0
                cx_ind2::Int64 = 1
                for i in opt[1].solver.Implicit_Options.gU_loc
                    rhs[cx_ind2] = sum(pmid[:].*c[i].cv_grad[:])+opt[1].solver.Implicit_Options.gU[i]-c[i].cv
                    cx_ind2 += 1
                end
                for i in opt[1].solver.Implicit_Options.gL_loc
                    rhs[cx_ind2] = sum(-pmid[:].*c[i].cc_grad[:])-opt[1].solver.Implicit_Options.gL[i]+c[i].cc
                    cx_ind2 += 1
                end
            end
            #println("ran to me 7:")
            #println("f.cv_grad: $(f.cv_grad)")
            #println("l: $(l)")
            #println("u: $(u)")
            #println("dcdx: $(dcdx)")
            model = buildlp([f.cv_grad[i] for i=1:np], dcdx, '<', rhs, l, u, opt[1].solver.LP_solver)
            result = solvelp(model)
            #println("result: $result.status")
            # Unpacks the results from the LP solver
            if (result.status == :Optimal)
                #println("ran me 8b")
                val::Float64 = result.objval + f_cv - sum([pmid[i]*f.cv_grad[i] for i=1:np])
                pnt::Vector{Float64} = vcat(mid.(Intv.(x_mc)),result.sol)
                feas = true
                mult::Vector{Float64} = result.attrs[:redcost]
                #println("val: $val")
                #println("pnt: $pnt")
                #println("feas: $feas")
                #println("mult: $mult")
                #println("ran me 8f")
            elseif (result.status == :Infeasible)
                #println("ran me 8b1")
                val = -Inf
                pnt = pmid
                feas = false
                mult = pmid
                #println("ran me 8f1")
            else
                #println("ran me 8b2")
                error("Solver error code $(result.status) in solver. Solution routine terminated.")
                #println("ran me 8b2")
            end

            mult_lo = [tol_eq(l[i],pnt[i],opt[1].solver.dual_tol) ? mult[i] : 0.0 for i=1:np]
            mult_hi = [tol_eq(u[i],pnt[i],opt[1].solver.dual_tol) ? mult[i] : 0.0 for i=1:np]
            #println("mult_lo: $mult_lo")
            #println("mult_hi: $mult_hi")
            temp = Any[mult_lo,mult_hi,val]
            #println("     -----     Lower LP Try Imp End     -----     ")
            return val, pnt, feas, temp
        catch
            # solve optimization problem via interval extension
            #println("ran catch!")
            FInt::MCInterval = opt[1].solver.Implicit_Options.f(Y[1:nx],Y[(nx+1):end])
            feas = true
            if (opt[1].solver.Implicit_Options.numConstr < 1)
            else
                GInt::Vector{MCInterval{Float64}} = opt[1].solver.Implicit_Options.g(Y[1:nx],Y[(nx+1):end])
                cInt::Vector{MCInterval{Float64}} = vcat(GInt[opt[1].solver.Implicit_Options.gU_loc]-opt[1].solver.Implicit_Options.gU[opt[1].solver.Implicit_Options.gU_loc],
                                                      -GInt[opt[1].solver.Implicit_Options.gL_loc]+opt[1].solver.Implicit_Options.gL[opt[1].solver.Implicit_Options.gL_loc])
                # Update To sparse later
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
            #println("     -----     Lower LP Catch Imp End     -----     ")
            return val, pnt, feas, temp
        end
end
