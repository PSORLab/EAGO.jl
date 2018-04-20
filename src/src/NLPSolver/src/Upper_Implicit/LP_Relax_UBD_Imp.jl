"""
    LP_Relax_UBD

Solves a upper bounding problem (LP) constructed by relaxing using nonsmooth
McCormick operators presented by McCormick1978, Mitsos2008. The relaxed problem
is then solved via an LP solvers. LP solver options and domain reduction
settings are controlled in the B&B algorithm using the global_options type.
Inputs:
* `X::Vector{Interval}`: Node over which to solve the lower problem
* `k::Int64`: Number of Iterations The B&B Algorithm has taken
* `pos::Int64`: Depth of Node in B&B Tree
* `opt::Any`: Option type containing problem information
* `temp::Any`: Information passed from lower problem to upper problem
Returns a tuple `(val,pnt,feas,X,[])` where
* `val::Float64`: Lower bound calculated
* `pnt::Vector{Float64}`: An array of length equal to X that gives the optimal
                          solution of the lower bound problem.
* `feas::Bool`: Returns true if the problem is feasible and false if it is infeasible
* `X::Vector{Interval}`: Node over which the problem was solved after potential
                         use of domain reduction techniques
"""
function LP_Relax_UBD_Imp(Y::Vector{Interval{Float64}},k::Int64,pos::Int64,opt::Any,temp::Any)

    nx::Int64 = opt[1].solver.Implicit_Options.nx
    np::Int64 = opt[1].numVar - nx

    l::Vector{Float64} = [Y[nx+i].lo for i=1:np]
    u::Vector{Float64} = [Y[nx+i].hi for i=1:np]
    pmid::Vector{Float64} = (l + u)/2.0

    param = InGenExpansionParams(opt[1].solver.Implicit_Options.h,
                                   opt[1].solver.Implicit_Options.hj,
                                   Y[1:nx],Y[(nx+1):opt[1].numVar],pmid,
                                   opt[1].solver.Implicit_Options.opts)
    x_mc = param[end]
    p_mc::Vector{SMCg{np,Float64}} = [SMCg{np,Float64}(pmid[i],
                                                       pmid[i],
                                                       seed_g(Float64,i,np),
                                                       seed_g(Float64,i,np),
                                                       Y[nx+i],
                                                       false,
                                                       Y[(nx+1):(nx+np)],
                                                       pmid) for i=1:np]
    # relaxation of function
    f::SMCg{np,Float64} = opt[1].solver.Implicit_Options.f(x_mc,p_mc)
    f_cc::Float64 = f.cc

    if opt[1].solver.Implicit_Options.numConstr>0
        c::Vector{SMCg{np,Float64}} = opt[1].solver.Implicit_Options.g(x_mc,p_mc)
        dcdx::SparseMatrixCSC{Float64,Int64} = spzeros(length(opt[1].solver.Implicit_Options.numConstr.gL_loc)+length(opt[1].solver.Implicit_Options.numConstr.gU_loc),np)
    else
        dcdx = spzeros(1,np)
    end

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

    # forms rhs of linear relaxations
    if opt[1].solver.Implicit_Options.numConstr>0
        rhs::Vector{Float64} = zeros(Float64,length(opt[1].solver.Implicit_Options.gL_loc)+length(opt[1].solver.Implicit_Options.gU_loc))
    else
        rhs = zeros(Float64,1)
    end

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

    model = buildlp([-f.cc_grad[i] for i=1:np], dcdx, '<', rhs, l, u, opt[1].solver.LP_solver)
    result = solvelp(model)

    if (result.status == :Optimal)
        val::Float64 = f_cc - result.objval - sum([pmid[i]*f.cc_grad[i] for i=1:np])
        pnt::Vector{Float64} = result.sol
        feas = true
    elseif (result.status == :Infeasible)
        val = -Inf
        pnt = pmid
        feas = false
    else
        error("Solver error code $(result.status) in $solver. Solution routine terminated.")
    end

    return val, pnt, feas, Any[feas,val]
end
