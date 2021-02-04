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

function set_tolerance_inner!(t::DefaultExt, alg, s, m::JuMP.Model, abs_tol::Float64)
    optimizer_name = MOI.SolverName()
    if optimizer_name === "EAGO: Easy Advanced Global Optimization"
        set_optimizer_attribute(m, "absolute_value", abs_tol)
    elseif optimizer_name === "SCIP"
        set_optimizer_attribute(m, "limits/absgap", abs_tol)
    elseif optimizer_name === "Alpine"
        set_optimizer_attribute(m, "absgap", abs_tol)
    elseif optimizer_name === "BARON"
        set_optimizer_attribute(m, "EpsA", abs_tol)
    elseif optimizer_name === "GAMS"
        set_optimizer_attribute(m, "OptCA", abs_tol)
    else
        error("A custom set_tolerance! function for solver = $optimizer_name
               specified for use with subproblem = $s in algorithm = $alg with
               extension = $t has not been defined and the selected solver is
               not support by default. Please open an issue requesting this feature
               at the following link https://github.com/PSORLab/EAGO.jl/issues.
               Extending the EAGO.set_tolerance! with a custom extension type
               will resolve this issue.")
    end
    return nothing
end

for (typ, fd) in SUBPROB_SYM
    @eval function load!(t::DefaultExt, alg::SIPResRev, s::$t, m::JuMP.Model, sr::SIPSubResult, i::Int)
        set_tolerance_inner!(t, alg, s, m, sr.$fd.tol, i)
        return nothing
    end
end

function sip_solve!(t, alg::SIPResRev, buffer::SIPSubResult, prob::SIPProblem,
                    result::SIPResult, cb::SIPCallback)

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
    print_summary!(LowerProblem(), verb, buffer)

    # solve inner program  and update lower discretization set
    is_llp1_nonpositive = true
    for i = 1:prob.nSIP
        sipRes_llp!(t, alg, LowerLevel1(), result, buffer, prob, cb, i)
        is_llp1_nonpositive &= buffer.objective_value > 0.0
        buffer.lbd_disc[i] .= buffer.pbar
        print_summary!(LowerLevel1(), verb, buffer, i)
    end
    # TODO: ADD TOLERANCE UPDATE CASE

    # if the lower problem is feasible then it's solution is the optimal value
    if is_llp1_nonpositive
        result.upper_bound = buffer.obj_value_lbd
        result.xsol .= buffer.lbd_x
        result.feasibility = true
        @goto main_end
    end
    push!(prob.lower_disc, deepcopy(buffer.lower_disc))

    # solve upper bounding problem, if feasible solve lower level problem,
    # and potentially update upper discretization set
    sipRes_bnd!(t, alg, UpperProblem(), buffer, eps_g, result, prob, cb)
    print_summary!(UpperProblem(), verb, buffer)
    if buffer.is_feasible_ubd
        is_llp2_nonpositive = true
        for i = 1:prob.nSIP
            sipRes_llp!(t, alg, LowerLevel2(), result, buffer, prob, cb, i)
            print_summary!(LowerLevel2(), verb, buffer, i)
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
        # TODO: ADD TOLERANCE UPDATE CASE
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
