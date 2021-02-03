# Copyright (c) 2018: Matthew Wilhelm & Matthew Stuber.
# This code is licensed under MIT license (see LICENSE.md for full details)
#############################################################################
# EAGO
# A development environment for robust and global optimization
# See https://github.com/PSORLab/EAGO.jl
#############################################################################
# src/eago_semiinfinite/algorithms/sip_hybrid.jl
# Defines the SIP-res algorithm which implements Algorithm #1 of XXX.
#############################################################################

struct SIPRes <: AbstractSIPAlgo end

function sip_solve!(t, alg::SIPRes, init_bnd, prob::SIPProblem, result::SIPResult, cb::SIPCallback)

    # initializes solution
    verb = prob.verbosity
    eps_g =  prob.initial_eps_g
    r = prob.initial_r

    @label main_iteration
    check_convergence(result, prob.absolute_tolerance, verb) && @goto main_end

    # solve lower bounding problem and check feasibility
    sipRes_bnd!(t, alg, LowerProblem(), buffer, 0.0, result, prob, cb)
    result.lower_bound = buffer.obj_value_lbd
    if buffer.is_feasible_lbd
        result.feasibility = false
        println("Terminated: lower bounding problem infeasible.")
        @goto main_end
    end
    print_summary!(verb, result, LowerProblem())

    # solve inner program  and update lower discretization set
    is_llp1_nonpositive = true
    for i = 1:prob.nSIP
        sipRes_llp!(t, alg, LowerLevel1(), result, buffer, prob, cb, i)
        is_llp1_nonpositive &= buffer.objective_value > 0.0
        buffer.lbd_disc[i] .= buffer.pbar
        ((verb == 1) || (verb == 2)) && print_summary!(buffer, "Lower LLP$i")
    end
    push!(prob.lower_disc, deepcopy(buffer.lower_disc))

    # if the lower problem is feasible then it's solution is the optimal value
    if is_llp1_nonpositive
        result.upper_bound = buffer.obj_value_lbd
        result.xsol .= buffer.lbd_x
        result.feasibility = true
        @goto main_end
    end

    # solve upper bounding problem, if feasible solve lower level problem,
    # and potentially update upper discretization set
    sipRes_bnd!(t, alg, UpperProblem(), buffer, eps_g, result, prob, cb)
    print_summary!(verb, buffer, :x, "UBD")
    if buffer.is_feasible_ubd
        is_llp2_nonpositive = true
        for i = 1:prob.nSIP
            sipRes_llp!(t, alg, LowerLevel2(), result, buffer, prob, cb, i)
            print_summary!(verb, buffer, :p, "Upper LLP$i")
            is_llp2_nonpositive &= buffer.objective_value > 0.0
            buffer.upper_disc[i] .= buffer.pbar
        end
        if is_llp2_nonpositive
            if buffer.obj_value_ubd <= result.upper_bound
                result.upper_bound = buffer.obj_value_ubd
                result.xsol .= buffer.ubd_x
            end
            eps_g /= r
        else
            push!(prob.upper_disc, deepcopy(buffer.upper_disc))
        end
    else
        eps_g /= r
    end

    # print iteration information and advance
    print_int!(verb, prob, k, result, eps_g, r)
    result.iteration_number += 1
    result.iteration_number < prob.iteration_limit && @goto main_iteration

    @label main_end
    return nothing
end
