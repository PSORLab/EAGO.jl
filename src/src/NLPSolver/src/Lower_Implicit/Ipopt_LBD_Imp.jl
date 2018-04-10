function LP_Relax_LBD_Imp(X::Vector{Interval{Float64}},
                          k::Int64,
                          pos::Int64,
                          opt::Any,
                          UBD::Float64)


        # feasibility-based bound tightening & Domain Reduction Presolves
        feas::Bool = true
        composite_DR_pre!(X,feas,k,pos,opt,UBD)

        # Solve the relaxed optimization problem
        if (feas == true)
            nx::Int64 = opt[1].solver.ImplicitOpts.nx
            np::Int64 = opt[1].numVar - nx
            p_L::Vector{Float64} = [X[i].lo for i=(nx+1):(nx+np)]
            p_U::Vector{Float64} = [X[i].hi for i=(nx+1):(nx+np)]
            p0::Vector{Float64} = (x_L + x_U)/2.0
            p_SMC::Vector{SMCg{np,Float64}} = [SMCg{np,Float64}(p0[i],
                                                                p0[i],
                                                                seed_g(np,i),
                                                                seed_g(np,i),
                                                                X[nx+i],
                                                                false,
                                                                X[(nx+1):(nx+np)],
                                                                p0) for i=1:np]
            # relaxation of function
            out,z,x_mc = GenExpansionParams(opt[1].solver.ImplicitOpts.h,
                                            opt[1].solver.ImplicitOpts.hj,
                                            X[1:nx],
                                            X[(nx+1):(nx+np)],
                                            p0,
                                            opt[1].solver.ImplicitOpts)

            f = opt[1].solver.ImplicitOpts.f(x_mc,p_SMC)
            f_cv::Float64 = f.cv

            if (opt[1].numConstr>0)

                c = opt[1].solver.ImplicitOpts.g(x_mc,p_SMC)

                c_cv::Vector{Float64} = [c[i].cv for i=1:opt[1].numConstr]
                c_cc::Vector{Float64} = [c[i].cc for i=1:opt[1].numConstr]
                dcdx_cv::Array{Float64,2} = Float64[c[i].cv_grad[j] for i=1:opt[1].numConstr,j=1:np]
                dcdx_cc::Array{Float64,2} = Float64[-c[i].cc_grad[j] for i=1:opt[1].numConstr,j=1:np]
                cval = vcat(c_cv[opt[1].gU_loc]-opt[1].gU[opt[1].gU_loc],-c_cc[opt[1].gL_loc]+opt[1].gL[opt[1].gL_loc])
                dcdx::Array{Float64,2} = vcat(dcdx_cv[opt[1].gU_loc,:],dcdx_cc[opt[1].gL_loc,:])

                rhs1::Vector{Float64} = [sum([x0[j]*dcdx_cv[i,j] for j=1:np]) for i=1:opt[1].numConstr]
                rhs2::Vector{Float64} = [sum([-x0[j]*dcdx_cv[i,j] for j=1:np]) for i=1:opt[1].numConstr]

                rhs = vcat(rhs1[opt[1].gU_loc]+opt[1].gU[opt[1].gU_loc]-c_cv[opt[1].gU_loc],
                           rhs2[opt[1].gL_loc]-opt[1].gL[opt[1].gL_loc]+c_cc[opt[1].gL_loc])
                # stores output for upper problem
                temp = [f,c]

                # sets up lp and solves
                temp_model = buildlp([f.cv_grad[i] for i=1:np], dcdx, '<', rhs, p_L, p_U, opt[1].solver.LP_solver)
                result = solvelp(temp_model)
            else

                # stores output for upper problem
                temp = [f_mc,[]]

                # sets up lp and solves
                temp_model = buildlp([f_mc.cv_grad[i] for i=1:np], zeros(np,np), '<', zeros(np), p_L, p_U, opt[1].solver.LP_solver)
                result = solvelp(temp_model)
            end

            # Unpacks the results from the LP solver
            if (result.status == :Optimal)
                val::Float64 = result.objval + f_cv - sum([x0[i]*f.cv_grad[i] for i=1:np])
                pnt::Vector{Float64} = result.sol
                feas = true
                mult::Vector{Float64} = result.attrs[:redcost]
            elseif (result.status == :Infeasible)
                val = -Inf
                pnt = x0
                feas = false
                mult = x0
            else
                error("Solver error code $(result.status) in $solver. Solution routine terminated.")
            end
        end

        # Formats Duality-Based Multipliers
        mult_lo = [tol_eq(X[i].lo,pnt[i],1E-12) ? mult[i] : 0.0 for i=1:np]
        mult_hi = [tol_eq(X[i].hi,pnt[i],1E-12) ? mult[i] : 0.0 for i=1:np]

        # Formats Runs Post Contractor
        composite_DR_post!(X,k,pos,opt[1],val,UBD,mult_lo,mult_hi)

    return val, pnt, feas, X, temp
end
