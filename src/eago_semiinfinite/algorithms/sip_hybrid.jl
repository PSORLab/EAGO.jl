# Copyright (c) 2018: Matthew Wilhelm & Matthew Stuber.
# This code is licensed under MIT license (see LICENSE.md for full details)
#############################################################################
# EAGO
# A development environment for robust and global optimization
# See https://github.com/PSORLab/EAGO.jl
#############################################################################
# src/eago_semiinfinite/algorithms/sip_hybrid.jl
# Defines the SIP-hybrid algorithm which implements Algorithm #2 of Djelassi,
# Hatim, and Alexander Mitsos. "A hybrid discretization algorithm with guaranteed
# feasibility for the global solution of semi-infinite programs."
# Journal of Global Optimization 68.2 (2017): 227-253.
#############################################################################

struct SIPHybrid <: AbstractSIPAlgo end

# Pseudocode currrently...
function sip_hybrid()

    @label main_algo
    # Solve LBD to get fLBD- at xLBD- and LBD <- fLBD-
    sipRes_bnd!(t, alg, LowerProblem(), buffer, 0.0, result, prob, cb)
    set_global_bound!(LowerProblem(), result, buffer)
    if !is_feasible(buffer)
        result.feasibility = false
        println("Terminated: lower bounding problem infeasible."); @goto main_end
    end
    print_summary!(verb, result, LowerProblem())

    # g+(xLBD-) and ybar <- Solve LLP with: xLBD-
    # solve inner program  and update lower discretization set
    is_llp1_nonpositive = true
    for i = 1:prob.nSIP
        sipRes_llp!(t, alg, LowerLevel1(), result, buffer, prob, cb, i)
        if buffer.objective_value > 0.0
            buffer.lower_disc[i] .= buffer.pbar
        ((verb == 1) || (verb == 2)) && print_summary!(buffer, "Lower LLP$i")
    end
    push!(prob.lower_disc, deepcopy(buffer.lower_disc))

    # if the lower problem is feasible then it's solution is the optimal value
    if is_llp1_nonpositive
        result.upper_bound = buffer.objective_value
        result.xsol .= buffer.xbar
        result.feasibility = true
        @goto main_end
    end
    if g(xLBD-, ybar) > 0.0
        # populate the discretization set with ybar
    else
        # update the parameter, ϵLLP
        ϵLLP = (g+(xLBD-) - g(xLBD-, ybar))/rLLP
    end

    @label upper_problem
    # Solve UBD to get fUBD- at xUBD-
    if is_feasible(#UBD)
        if g+(xUBD-) <= 0.0
            if fUBD- < UBD
                UDB <- fUBD
                xstar <- xUBD-
            end
                ϵg /= rg
                break
            else
                # populate the discretization set with ybar
            end
        else
            ϵg /= rg
        end
    end

    # check for termination...
    l = 0
    @label oracle_loop
        # SOLVE RES
        if η+ < 0
            LBD <= fRES
            fRES <= 0.5*(LBD + UBD)
            continue
        elseif ηbar > 0
            if
            elseif
            else
                @goto main_iteration
            end
        else
            @goto main_iteration
        end
    end

    return result
end
