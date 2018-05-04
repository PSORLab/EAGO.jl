function NLPSolve_UBD(X,k::Int64,pos::Int64,opts,temp)

    # reloads upper model
    mod = deepcopy(opts[1].UBDmodel)
    loadproblem!(mod,opts[1].numVar,opts[1].numConstr,
                [X[i].lo for i=1:length(X)],[X[i].hi for i=1:length(X)],
                opts[1].gL,opts[1].gU,opts[1].sense,opts[1].d)

    #setwarmstart!(opts[1].UBDmodel,mid.(X))

    # solves upper model
    optimize!(mod)
    out = status(mod)

    pnt::Vector{Float64} = getsolution(mod)
    val::Float64 = getobjval(mod)
    if (out == :Optimal)
        feas::Bool = true
    elseif (out == :Infeasible)
        feas = false
    else
        val = (opts[1].f(X)).hi
        pnt = mid.(X)
        feas = true
        #error("Solver error code $status in upper problem. Solution routine terminated.")
    end


    return val, pnt, feas, Any[feas,val]
end
