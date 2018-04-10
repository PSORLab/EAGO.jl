function NLPSolve_UBD(X,k::Int64,pos::Int64,opts,temp)

    # reloads upper model
    loadproblem!(opts[1].UBDmodel,opts[1].numVar,opts[1].numConstr,
                [X[i].lo for i=1:length(X)],[X[i].hi for i=1:length(X)],
                opts[1].gL,opts[1].gU,opts[1].sense,opts[1].d)

    #setwarmstart!(opts[1].UBDmodel,mid.(X))

    # solves upper model
    optimize!(opts[1].UBDmodel)
    out = status(opts[1].UBDmodel)

    pnt::Vector{Float64} = getsolution(opts[1].UBDmodel)
    val::Float64 = getobjval(opts[1].UBDmodel)
    if (out == :Optimal)
        println("feas true")
        feas::Bool = true
    elseif (out == :Infeasible)
        feas == false
        println("feas false")
    else
        val = (opts[1].f(X)).hi
        pnt = mid.(X)
        feas = true
        #error("Solver error code $status in upper problem. Solution routine terminated.")
    end


    return val, pnt, feas, Any[feas,val]
end
