# Copyright (c) 2018: Matthew Wilhelm, Robert Gottlieb, Dimitri Alston,
# Matthew Stuber, and the University of Connecticut (UConn).
# This code is licensed under the MIT license (see LICENSE.md for full details).
################################################################################
# EAGO
# A development environment for robust and global optimization.
# https://github.com/PSORLab/EAGO.jl
################################################################################
# src/eago_semiinfinite/nonconvex_algorithms/sip_hybrid.jl
# Defines the SIP-res algorithm which implements Algorithm #1 of XXX.           #TODO: Find the citation for this.
################################################################################

"""
    SIPRes

Specifies that the SIPRes algorithm which implements Algorithm #1 of Djelassi,
Hatim, and Alexander Mitsos. "A hybrid discretization algorithm with guaranteed
feasibility for the global solution of semi-infinite programs." Journal of Global
Optimization 68.2 (2017): 227-253 should be used.
"""
struct SIPRes <: AbstractSIPAlgo end

get_sip_optimizer(t::DefaultExt, alg::SIPRes, s::S) where S <: AbstractSubproblemType = EAGO.Optimizer

function get_disc_set(t::ExtensionType, alg::SIPRes, s::LowerProblem, sr::SIPSubResult, i::Int)
    sr.disc_l[i]
end
function get_disc_set(t::ExtensionType, alg::SIPRes, s::UpperProblem, sr::SIPSubResult, i::Int)
    sr.disc_u[i]
end

function sip_solve!(t::ExtensionType, alg::SIPRes, buffer::SIPSubResult, prob::SIPProblem,
                    result::SIPResult, cb::SIPCallback)

    verb = prob.verbosity

    # Initializes solution
    @label main_iteration
    check_convergence(result, prob.abs_tolerance, verb) && @goto main_end

    # Solve lower bounding problem and check feasibility
    sip_bnd!(t, alg, LowerProblem(), buffer, result, prob, cb)
    result.lower_bound = buffer.lbd.obj_val
    if !buffer.lbd.feas
        result.feasibility = false
        println("Terminated: lower bounding problem infeasible.")
        @goto main_end
    end
    print_summary!(LowerProblem(), verb, buffer)

    # Solve inner program and update lower discretization set
    is_llp1_nonpositive = true
    for i = 1:prob.nSIP
        sip_llp!(t, alg, LowerLevel1(), result, buffer, prob, cb, i)
        buffer.disc_l_buffer .= buffer.llp1.sol
        print_summary!(LowerLevel1(), verb, buffer, i)
        if buffer.llp1.obj_val <= 0.0
            continue
        else
            push!(buffer.disc_l[i], deepcopy(buffer.disc_l_buffer))
            is_llp1_nonpositive = false
        end
    end

    # If the lower problem is feasible, then it's solution is the optimal value
    if is_llp1_nonpositive
        result.upper_bound = buffer.lbd.obj_val
        result.xsol .= buffer.lbd.sol
        result.feasibility = true
        @goto main_end
    end

    # Solve upper bounding problem, if feasible solve lower level problem,
    # and potentially update upper discretization set
    sip_bnd!(t, alg, UpperProblem(), buffer, result, prob, cb)
    print_summary!(UpperProblem(), verb, buffer)
    if buffer.ubd.feas
        is_llp2_nonpositive = true
        for i = 1:prob.nSIP
            sip_llp!(t, alg, LowerLevel2(), result, buffer, prob, cb, i)
            buffer.disc_u_buffer .= buffer.llp2.sol
            print_summary!(LowerLevel2(), verb, buffer, i)
            if buffer.llp2.obj_val <= 0.0
                buffer.eps_g[i] /= buffer.r_g
                continue
            else
                push!(buffer.disc_u[i], deepcopy(buffer.disc_u_buffer))
                is_llp2_nonpositive = false
            end
        end
        if is_llp2_nonpositive
            if buffer.ubd.obj_val <= result.upper_bound
                result.upper_bound = buffer.ubd.obj_val
                result.xsol .= buffer.ubd.sol
            end
        end
    else
        buffer.eps_g ./= buffer.r_g
    end

    println("FINISHED ONE ITERATION")

    # Print iteration information and advance
    print_int!(verb, prob, result, buffer.r_g)
    result.iteration_number += 1
    result.iteration_number < prob.iteration_limit && @goto main_iteration

    @label main_end
    return nothing
end
