# Copyright (c) 2018: Matthew Wilhelm & Matthew Stuber.
# This code is licensed under MIT license (see LICENSE.md for full details)
#############################################################################
# EAGO
# A development environment for robust and global optimization
# See https://github.com/PSORLab/EAGO.jl
#############################################################################
# src/eago_semiinfinite/algorithms/sip_hybrid.jl
# Defines the revised SIP-res algorithm which implements Algorithm #1 of Djelassi,
# Hatim, and Alexander Mitsos. "A hybrid discretization algorithm with guaranteed
# feasibility for the global solution of semi-infinite programs."
# Journal of Global Optimization 68.2 (2017): 227-253.
#############################################################################

struct SIPResRev <: AbstractSIPAlgo end

function load!(t::DefaultExt, alg::SIPResRev, s::LowerLevel1, m::JuMP.Model, sr::SIPSubResult, i::Int)
    set_tolerance_inner!(t, alg, s, m, sr.eps_l[i])
end
function load!(t::DefaultExt, alg::SIPResRev, s::LowerLevel2, m::JuMP.Model, sr::SIPSubResult, i::Int)
    set_tolerance_inner!(t, alg, s, m, sr.eps_u[i])
end
function load!(t::DefaultExt, alg::SIPResRev, s::LowerProblem, m::JuMP.Model, sr::SIPSubResult, i::Int)
    set_tolerance_inner!(t, alg, s, m, sr.lbd.tol)
end
function load!(t::DefaultExt, alg::SIPResRev, s::UpperProblem, m::JuMP.Model, sr::SIPSubResult, i::Int)
    set_tolerance_inner!(t, alg, s, m, sr.ubd.tol)
end

function get_disc_set(t::ExtensionType, alg::SIPResRev, s::LowerProblem, sr::SIPSubResult, i::Int)
    sr.disc_l[i]
end
function get_disc_set(t::ExtensionType, alg::SIPResRev, s::UpperProblem, sr::SIPSubResult, i::Int)
    sr.disc_u[i]
end

function sip_solve!(t, alg::SIPResRev, buffer::SIPSubResult, prob::SIPProblem,
                    result::SIPResult, cb::SIPCallback)

    verb = prob.verbosity

    # initializes solution
    @label main_iteration
    check_convergence(result, prob.absolute_tolerance, verb) && @goto main_end

    # solve lower bounding problem and check feasibility
    sip_bnd!(t, alg, LowerProblem(), buffer, result, prob, cb)
    result.lower_bound = buffer.obj_value_lbd
    if buffer.is_feasible_lbd
        result.feasibility = false
        println("Terminated: lower bounding problem infeasible.")
        @goto main_end
    end
    print_summary!(LowerProblem(), verb, buffer)

    # solve inner program and update lower discretization set
    is_llp1_nonpositive = true
    for i = 1:prob.nSIP
        sip_llp!(t, alg, LowerLevel1(), result, buffer, prob, cb, i)
        buffer.disc_l_buffer .= buffer.pbar
        print_summary!(LowerLevel1(), verb, buffer, i)
        if buffer.llp1.obj_bnd <= 0.0
            continue
        elseif buffer.llp1.obj_val > 0.0
            push!(prob.disc_l[i], deepcopy(buffer.disc_l_buffer))
            is_llp1_nonpositive = false
        else
            buffer.eps_l[i] = (buffer.llp1.obj_bnd - buffer.llp1.obj_val)/buffer.r_l
            is_llp1_nonpositive = false
        end
    end

    # if the lower problem is feasible then it's solution is the optimal value
    if is_llp1_nonpositive
        result.upper_bound = buffer.lbd.obj_val
        result.xsol .= buffer.lbd.sol
        result.feasibility = true
        @goto main_end
    end

    # solve upper bounding problem, if feasible solve lower level problem,
    # and potentially update upper discretization set
    sip_bnd!(t, alg, UpperProblem(), buffer, result, prob, cb)
    print_summary!(UpperProblem(), verb, buffer)
    if buffer.is_feasible_ubd
        is_llp2_nonpositive = true
        for i = 1:prob.nSIP
            sip_llp!(t, alg, LowerLevel2(), result, buffer, prob, cb, i)
            buffer.disc_u_buffer[i] .= buffer.pbar
            print_summary!(LowerLevel2(), verb, buffer, i)
            if buffer.llp2.obj_bnd <= 0.0
                buffer.eps_g[i] /= buffer.r_g
                continue
            else
                push!(prob.disc_u[i], deepcopy(buffer.disc_u_buffer))
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

    # print iteration information and advance
    print_int!(verb, prob, k, result, buffer.r_g)
    result.iteration_number += 1
    result.iteration_number < prob.iteration_limit && @goto main_iteration

    @label main_end
    return nothing
end
