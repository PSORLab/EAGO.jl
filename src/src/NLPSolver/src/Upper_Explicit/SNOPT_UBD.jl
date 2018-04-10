"""
    snopt_callback_UBD

Callback function used by SNOPT in the upper bounding problem (local NLP). Inputs
are given by:
* `y::Vector{Float64}` Point to evaluate in X
* `opt` - Option type containing problem information
Returns the tuple `(f, c, dfdx, dcdx, fail)` where
* `f::Float64`: Value of the objective function
* `c::Array{Float64,1}`: Value of the constraint function
* `dfdx::Array{Float64,1}`: Gradient of the objective function
* `dcdx::Array{Float64,2}`: Jacobian of the constraint function
* `fail::Bool`: Flag used by SNOPT solver to indicate if SNOPT solver failed.
"""
function snopt_callback_UBD(y::Vector{Float64},opt::EAGO_Inner_NLP)

    f::Float64 = opt.f(y)
    dfdx::Vector{Float64} = ForwardDiff.gradient(opt.f,y)
    if (opt.numConstr>0)
        g_eval = opt.g(y)
        j_eval = ForwardDiff.jacobian(x->opt.g(x),y)
        c::Vector{Float64} = vcat(g_eval[opt.gU_loc]-opt.gU[opt.gU_loc],-g_eval[opt.gL_loc]+opt.gL[opt.gL_loc])
        dcdx::Array{Float64,2} = vcat(j_eval[opt.gU_loc,:],-j_eval[opt.gL_loc,:])
    else
        c = [-1.0]
        dcdx = zeros(1,opt.numVar)
    end
    fail = false

    return f, c, dfdx, dcdx, fail
end

"""
    SNOPT_UBD

Solves a upper bounding problem (local NLP) using SNOPT. SNOPT options and
domain reduction settings are controlled in the B&B algorithm using the
global_options type. Inputs are:
* `X::Vector{Interval{Float64}}`: Node over which to solve the lower problem
* `k::Int64`: Number of Iterations The B&B Algorithm has taken
* `pos::Int64`: Depth of Node in B&B Tree
* `opt`: Option type containing problem information
* `temp`: Info passed from lower problem to upper problem evaluation
Returns the tuple `(val,pnt,feas,X,Any[feas,val])` where
* `val::Float64`: Upper bound calculated
* `pnt::Array{Float64,1}`: An array of length equal to X that gives the
                           optimal solution of the upper bound problem.
* `feas::Bool`: Returns true if the problem is feasible and false if it is
                infeasible
* `Any[feas,val]`: Info for selective branching
"""
function SNOPT_UBD(X::Vector{Interval{Float64}},
                   k::Int64,
                   pos::Int64,
                   opt,
                   temp)
        # sets up problem
        feas::Bool = false
        x_Li::Vector{Float64} = [X[i].lo for i=1:opt[1].numVar]
        x_Ui::Vector{Float64} = [X[i].hi for i=1:opt[1].numVar]
        x0i::Vector{Float64} = (x_Li + x_Ui)/2.0
        options = Dict{String, Any}()
        options["Derivative option"] = 1
        options["Verify level"] = 0

        # solve problem and unpacks variables
        if (opt[1].solver.UBD_full_depth < pos)
            options["Major optimality tolerance"] = Inf
        else
            options["Major optimality tolerance"] = 1e-6
        end
        options["Major Feasibility Constraint"]
        TT = STDOUT
        redirect_stdout()
        pnt, val, status, mult = Snopt.snopt(y::Vector{Float64} -> opt[2].fg_SNOPT_UBD(y), x0i, x_Li, x_Ui, options)
        redirect_stdout(TT)

        if (status == 1 || status == 2)
            feas = true
        elseif (status == 11 || status == 12 || status == 13 || status == 14)
            println("status code: $status")
            feas == false
        else
            error("Solver error code $status in Snopt. Solution routine terminated.")
        end
        # output
        return val, pnt, feas, Any[feas,val]
end

function SNOPT_UBD(X::Vector{MCInterval{Float64}},
                   k::Int64,
                   pos::Int64,
                   opt,
                   temp)
        # sets up problem
        feas::Bool = false
        x_Li::Vector{Float64} = [X[i].lo for i=1:opt[1].numVar]
        x_Ui::Vector{Float64} = [X[i].hi for i=1:opt[1].numVar]
        x0i::Vector{Float64} = (x_Li + x_Ui)/2.0
        options = Dict{String, Any}()
        options["Derivative option"] = 1
        options["Verify level"] = 0

        # solve problem and unpacks variables
        if (opt[1].solver.UBD_full_depth < pos)
            options["Major optimality tolerance"] = Inf
        else
            options["Major optimality tolerance"] = 1e-6
        end
        TT = STDOUT
        redirect_stdout()
        pnt, val, status, mult = Snopt.snopt(y::Vector{Float64} -> opt[2].fg_SNOPT_UBD(y), x0i, x_Li, x_Ui, options)
        redirect_stdout(TT)

        if (status == 1 || status == 2)
            feas = true
        elseif (status == 11 || status == 12 || status == 13 || status == 14)
            println("status code: $status")
            feas == false
        else
            error("Solver error code $status in Snopt. Solution routine terminated.")
        end
        # output
        return val, pnt, feas, Any[feas,val]
end
