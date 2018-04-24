#=
function Quadratic_ObjCV_Callback(x::Vector{Float64},X::Vector{Interval{Float64}},opt,a::Vector{Float64})
            x_SMC::Vector{SMC} = [SMC(x[i],x[i],X[i]) for i=1:opt.numVar]
            y::SMC = opt.f(x_SMC-a)
            return y.cv
end

function Quadratic_LBD(X::Vector{Interval{Float64}},k::Int64,pos::Int64,opt,UBD::Float64)

            # feasibility-based bound tightening & Domain Reduction Presolves
            feas::Bool = true
            val::Float64 = 0.0
            pnt::Vector{Float64} = zeros(Float64,opt[1].numVar)
            mult::Vector{Float64} = zeros(Float64,opt[1].numVar)
            composite_DR_pre!(X,feas,k,pos,opt,UBD)

            # takes hessian of shifted problem
            if (feas == true)
                        hess_r = ReverseDiff.hessian(x->opt[2].Quadratic_ObjCV_Callback(x,X,mid.(X)),zeros(Float64,opt[1].numVar))
                        grad_r = ReverseDiff.gradient(x->opt[2].Quadratic_ObjCV_Callback(x,X,mid.(X)),zeros(Float64,opt[1].numVar))
                        x_L::Vector{Float64} = [X[i].lo-mid(X[i]) for i=1:opt[1].numVar]
                        x_U::Vector{Float64} = [X[i].hi-mid(X[i]) for i=1:opt[1].numVar]
                        x0::Vector{Float64} = (x_L + x_U)/2.0
                        Xt::Vector{Float64} = [Interval(x_L,x_U) for i=1:opt[1].numVar]
                        x_SMC::Vector{SMCg{opt[1].numVar,Float64}} = [SMCg{opt[1].numVar,Float64}(x0[i],x0[i],seed_g(opt[1].numVar,i),seed_g(opt[1].numVar,i),X[i],false,Xt,x0) for i=1:opt[1].numVar]

                        f::SMCg{opt[1].numVar,Float64} = opt[1].f(x_SMC+mid(X[i]))
                        f_cv::Float64 = f.cv
                        if (opt[1].numConstr>0)
                           # relaxation of constraints
                           c::Vector{SMCg{opt[1].numVar,Float64}} = opt[1].g(x_SMC+mid.(X))
                           c_cv::Vector{Float64} = [c[i].cv for i=1:opt[1].numConstr]
                           dcdx_cv::Array{Float64,2} = [c[i].cv_grad[j] for i=1:opt[1].numConstr, j=1:opt[1].numVar]
                           rhs::Vector{Float64} = [sum([x0[j]*dcdx_cv[i,j] for j=1:opt[1].numVar]) for i=1:opt[1].numConstr]-c_cv
                           # stores output for upper problem
                           temp = f,c
                        else
                           # stores output for upper problem
                           temp = f,[]
                        end
                        result = quadprog(grad_r,hess_r,dcdx_cv,'<',rhs,x_L,x_U,IpoptSolver())

                        # Unpacks the results from the LP solver
                        if (result.status == :Optimal)
                           val = result.objval + f_cv
                           pnt = result.sol + mid.(X)
                           feas = true
                           mult = result.attrs[:redcost]
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
            mult_lo = [tol_eq(X[i].lo,pnt[i],1E-12) ? mult[i] : 0.0 for i=1:opt[1].numVar]
            mult_hi = [tol_eq(X[i].hi,pnt[i],1E-12) ? mult[i] : 0.0 for i=1:opt[1].numVar]

            # Formats Runs Post Contractor
            composite_DR_post!(X,k,pos,opt[1],val,UBD,mult_lo,mult_hi)

            return val, pnt, feas, X, temp
end
=#
