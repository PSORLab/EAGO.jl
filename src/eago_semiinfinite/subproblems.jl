# Copyright (c) 2018: Matthew Wilhelm, Robert Gottlieb, Dimitri Alston,
# Matthew Stuber, and the University of Connecticut (UConn).
# This code is licensed under the MIT license (see LICENSE.md for full details).
################################################################################
# EAGO
# A development environment for robust and global optimization.
# https://github.com/PSORLab/EAGO.jl
################################################################################
# src/eago_semiinfinite/subproblems.jl
# Defines utilities for generic SIP subroutines.
################################################################################

"""
    build_model

Create the model and variables used with extension `t::EAGO.ExtensionType` in algorithm
`a::AbstractSIPAlgo` in subproblem `s::AbstractSubproblemType` via the command
`build_model(t::ExtensionType, a::AbstractSIPAlgo, s::AbstractSubproblemType, p::SIPProblem)`.
"""
function build_model(t::DefaultExt, a::A, s::S, p::SIPProblem) where {A <: AbstractSIPAlgo, S <: AbstractSubproblemType}
    model = Model(get_sip_optimizer(t,a,s))
    for (k,v) in get_sip_kwargs(s,p)
        MOI.set(model, MOI.RawParameter(String(k)), v)
    end
    vL, vU, nv = get_bnds(s,p)
    @variable(model, vL[i] <= v[i=1:nv] <= vU[i])
    return model, v
end
function build_model(t::ExtensionType, a::A, s::S, p::SIPProblem) where {A <: AbstractSIPAlgo, S <: AbstractSubproblemType}
    build_model(DefaultExt(), a, s, p)
end

###
### General set tolerances
###
function set_tolerance_inner!(t::DefaultExt, alg, s, m::JuMP.Model, abs_tol::Float64)
    optimizer_name = JuMP.solver_name(m)
    if optimizer_name === "EAGO - Easy Advanced Global Optimization"
        set_optimizer_attribute(m, "absolute_tolerance", abs_tol)
        #set_optimizer_attribute(m, "constr_viol_tol", c_tol)
        #set_optimizer_attribute(m, "relative_tolerance", rel_tol)
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
function set_tolerance!(t::DefaultExt, alg::AbstractSIPAlgo, s::S, m::JuMP.Model,
                        sr::SIPSubResult, i::Int) where {S <: AbstractSubproblemType}
    return nothing
end
function set_tolerance!(t::ExtensionType, alg::A, s::S, m::JuMP.Model,
                        sr::SIPSubResult, i::Int) where {A <: AbstractSIPAlgo,
                                                         S <: AbstractSubproblemType}
    set_tolerance!(t, alg, s, m, sr, i)
end

###
### General get discretation set
###
function get_disc_set(t::DefaultExt, alg::AbstractSIPAlgo, s::S, sr::SIPProblem, i::Int) where {S <: AbstractSubproblemType}
    return Vector{Float64}[]
end
function get_disc_set(t::ExtensionType, alg::AbstractSIPAlgo, s::S, sr::SIPProblem, i::Int) where {S <: AbstractSubproblemType}
    get_disc_set(DefaultExt(), alg, s, p, i)
end

###
### Add constraints solely on uncertainty
###
function add_uncertainty_constraint!(m::JuMP.Model, prob::SIPProblem)
    #if !isnothing(model.polyhedral_uncertainty_set)
        #@constraint(model, )
    #end
    #if !isnothing(model.ellipsodial_uncertainty_set)
        #@constraint(model, )
    #end
    return nothing
end

get_xbar(t::DefaultExt, alg::AbstractSIPAlgo, s::LowerLevel1, sr::SIPSubResult) = sr.lbd.sol
get_xbar(t::DefaultExt, alg::AbstractSIPAlgo, s::LowerLevel2, sr::SIPSubResult) = sr.ubd.sol
get_xbar(t::DefaultExt, alg::AbstractSIPAlgo, s::LowerLevel3, sr::SIPSubResult) = sr.lbd.res

function llp_check(islocal::Bool, t::MOI.TerminationStatusCode, r::MOI.ResultStatusCode)
    flag = true
    if islocal && ((t != MOI.LOCALLY_SOLVED) && (t != MOI.ALMOST_LOCALLY_SOLVED))
        error("Lower problem did not solve to local optimality.")
    elseif t !== MOI.OPTIMAL
        error("Error in lower level problem. Termination status = $t, primal status = $r.")
    end
    return flag
end

"""
    sip_llp!

Solves the lower level problem for the ith-SIP used with extension `t::EAGO.ExtensionType`
in algorithm `a::AbstractSIPAlgo` in subproblem `s::AbstractSubproblemType` via
the command  `sip_llp!(t::ExtensionType, a::AbstractSIPAlgo, s::AbstractSubproblemType, ..., i, tol)`.
"""
function sip_llp!(t::DefaultExt, alg::A, s::S, result::SIPResult,
                  sr::SIPSubResult, prob::SIPProblem, cb::SIPCallback,
                  i::Int64, tol::Float64 = -Inf) where {A <: AbstractSIPAlgo, S <: AbstractSubproblemType}

    # Build the model
    m, p = build_model(t, alg, s, prob)
    set_tolerance!(t, alg, s, m, sr, i)

    # Define the objective
    xbar = get_xbar(t, alg, s, sr)
    g(p...) = cb.gSIP[i](xbar, p)
    register(m, :g, prob.np, g, autodiff=true)
    if isone(prob.np)
        nl_obj = :(g($(p[1])))
    else
        nl_obj = Expr(:call)
        push!(nl_obj.args, :g)
        for i in 1:prob.np
            push!(nl_obj.args, p[i])
        end
    end
    set_nonlinear_objective(m, MOI.MAX_SENSE, nl_obj)

    # Add uncertainty constraints
    add_uncertainty_constraint!(m, prob)

    # Optimize model and check status
    JuMP.optimize!(m)
    tstatus = JuMP.termination_status(m)
    rstatus = JuMP.primal_status(m)
    feas = llp_check(prob.local_solver, tstatus, rstatus)

    # Fill buffer with subproblem result info
    psol = JuMP.value.(p)
    load!(s, sr, feas, JuMP.objective_value(m), JuMP.objective_bound(m), psol)
    result.solution_time += MOI.get(m, MOI.SolveTimeSec())

    return nothing
end
function sip_llp!(t::ExtensionType, alg::A, s::S, result::SIPResult,
                     sr::SIPSubResult, prob::SIPProblem, cb::SIPCallback,
                     i::Int64) where {A <: AbstractSIPAlgo, S <: AbstractSubproblemType}
    sip_llp!(DefaultSubproblem(), alg, s, result, sr, prob, cb, i)
end

"""
    sip_bnd!

Solves the bounding problem for the ith-SIP used with extension `t::EAGO.ExtensionType`
in algorithm `a::AbstractSIPAlgo` in subproblem `s::AbstractSubproblemType` via
the command `sip_bnd!(t::ExtensionType, a::AbstractSIPAlgo, s::AbstractSubproblemType, ..., i, tol)`.
"""
function sip_bnd!(t::ExtensionType, alg::A, s::S, sr::SIPSubResult, result::SIPResult,
                  prob::SIPProblem, cb::SIPCallback) where {A <: AbstractSIPAlgo,
                                                            S <: AbstractSubproblemType}
    sip_bnd!(DefaultExt(), alg, s, sr, result, prob, cb)
end
function sip_bnd!(t::DefaultExt, alg::A, s::S, sr::SIPSubResult, result::SIPResult,
                  prob::SIPProblem, cb::SIPCallback) where {A <: AbstractSIPAlgo,
                                                            S <: AbstractSubproblemType}

    # Create JuMP model
    m, x = build_model(t, alg, s, prob)

    for i = 1:prob.nSIP
        ϵ_g = get_eps(s, sr, i)
        disc_set = get_disc_set(t, alg, s, sr, i)
        for j = 1:length(disc_set)
            gi = Symbol("g$i$j")
            g(x...) = cb.gSIP[i](x, disc_set[j])
            register(m, gi, prob.nx, g, autodiff=true)
            gic = Expr(:call)
            push!(gic.args, gi)
            for i in 1:prob.nx
                push!(gic.args, x[i])
            end
            JuMP.add_nonlinear_constraint(m, :($gic + $ϵ_g <= 0))
        end
    end

    # Define the objective
    obj(x...) = cb.f(x)
    register(m, :obj, prob.nx, obj, autodiff=true)
    if isone(prob.nx)
        nl_obj = :(obj($(x[1])))
    else
        nl_obj = Expr(:call)
        push!(nl_obj.args, :obj)
        for i in 1:prob.nx
            push!(nl_obj.args, x[i])
        end
    end
    set_nonlinear_objective(m, MOI.MIN_SENSE, nl_obj)

    # Optimize model and check status
    JuMP.optimize!(m)
    t_status = JuMP.termination_status(m)
    r_status = JuMP.primal_status(m)
    feas = bnd_check(prob.local_solver, t_status, r_status)

    # Fill buffer with subproblem result info
    load!(s, sr, feas, JuMP.objective_value(m), JuMP.objective_bound(m), JuMP.value.(x))
    result.solution_time += MOI.get(m, MOI.SolveTimeSec())

    return nothing
end

"""
    sip_res!

Solves the restriction problem for extension `t::EAGO.ExtensionType`
in algorithm `a::AbstractSIPAlgo` in subproblem `s::AbstractSubproblemType` via
the command `sip_res!(t::ExtensionType, a::AbstractSIPAlgo, ...)`.
"""
function sip_res!(t::ExtensionType, alg::A, sr::SIPSubResult, result::SIPResult,
                  prob::SIPProblem, cb::SIPCallback) where {A <: AbstractSIPAlgo}
    sip_res!(DefaultExt(), alg, sr, result, prob, cb)
end
function sip_res!(t::DefaultExt, alg::A, sr::SIPSubResult, result::SIPResult,
                  prob::SIPProblem, cb::SIPCallback) where {A <: AbstractSIPAlgo}

    # Create JuMP model and variables
    s = ResProblem()
    m, x = build_model(t, alg, s, prob)
    @variable(m, η)

    # Add discretized semi-infinite constraint
    for i = 1:prob.nSIP
        disc_set = get_disc_set(t, alg, s, sr, i)
        for j = 1:length(disc_set)
            gi = Symbol("g$i$j")
            g(x...) = cb.gSIP[i](x, disc_set[j])
            register(m, gi, prob.nx, g, autodiff=true)
            gic = Expr(:call)
            push!(gic.args, gi)
            for i in 1:prob.nx
                push!(gic.args, x[i])
            end
            gip = Expr(:call)
            push!(gip.args, :+)
            push!(gip.args, gic)
            push!(gip.args, η)
            JuMP.add_nonlinear_constraint(m, :($gip <= 0))
        end
    end

    # Add epigraph reformulated objective
    obj(x...) = cb.f(x)
    register(m, :f, prob.nx, obj, autodiff=true)
    if isfinite(sr.fRes)
        if isone(prob.nx)
            nl_obj = :(f($(x[1])))
        else
            nl_obj = Expr(:call)
            push!(nl_obj.args, :f)
            for i in 1:prob.nx
                push!(nl_obj.args, x[i])
            end
        end
        JuMP.add_nonlinear_constraint(m, :($nl_obj + $(sr.fRes) <= 0))
    end

    # Define the objective
    @objective(m, Min, -η)

    # Optimize model and check status
    JuMP.optimize!(m)
    t_status = JuMP.termination_status(m)
    r_status = JuMP.primal_status(m)
    feas = llp_check(prob.local_solver, t_status, r_status)

    # Fill buffer with subproblem result info
    load!(s, sr, feas, JuMP.objective_value(m), JuMP.objective_bound(m), JuMP.value.(x))
    result.solution_time += MOI.get(m, MOI.SolveTimeSec())

    return nothing
end

# Get optimizer for use in subproblems
function get_sip_optimizer(t::DefaultExt, alg::A, s::S) where {A<:AbstractSIPAlgo, S<:AbstractSubproblemType}
    return EAGO.Optimizer
end

"""
    get_sip_optimizer

Specifices the optimizer to be used in extension `t::EAGO.ExtensionType` with
algorithm `alg::AbstractSIPAlgo` in subproblem `s::AbstractSubproblemType` via the
command `get_sip_optimizer(t::ExtensionType, alg::AbstractSIPAlgo, s::AbstractSubproblemType)`.
"""
function get_sip_optimizer(t::ExtensionType, alg::A, s::AbstractSubproblemType) where {A <: AbstractSIPAlgo}
    return get_sip_optimizer(DefaultExt(), alg, s)
end

# Printing
function print_int!(verb::Int, prob::SIPProblem, result::SIPResult, r::Float64)

    k = result.iteration_number
    lbd = result.lower_bound
    ubd = result.upper_bound

    if (prob.verbosity == 1 || prob.verbosity == 2)

        # Prints header line every hdr_intv times
        if (mod(k, prob.header_interval) == 0 || k == 1)
            println("| Iteration | Lower Bound | Upper Bound |   r   |  Gap  |  Ratio  |")
        end

        # Prints iteration summary every prnt_intv times
        if mod(k, prob.print_interval) == 0

            print_str = "| "

            max_len = 15
            temp_str = string(k)
            len_str = length(temp_str)
            print_str *= (" "^(max_len - len_str))*temp_str*" | "

            max_len = 15
            lower_adj = lbd
            temp_str = string(round(lower_adj, digits = 6))
            len_str = length(temp_str)
            print_str *= (" "^(max_len - len_str))*temp_str*" | "

            max_len = 15
            upper_adj = ubd
            temp_str = string(round(upper_adj, digits = 6))
            len_str = length(temp_str)
            print_str *= (" "^(max_len - len_str))*temp_str*" | "

            max_len =15
            temp_str = string(round(r, digits = 5))
            len_str = length(temp_str)
            print_str *= (" "^(max_len - len_str))*temp_str*" | "

            max_len = 15
            temp_str = string(round((ubd-lbd), digits = 5))
            len_str = length(temp_str)
            print_str *= (" "^(max_len - len_str))*temp_str*" | "

            max_len = 18
            temp_str = string(round((ubd-lbd)/abs(ubd), digits = 6))
            len_str = length(temp_str)
            print_str *= (" "^(max_len - len_str))*temp_str*" |"

            println(print_str)
        end
    end
    return nothing
end

function summary_inner!(verb::Int64, val::Float64, x::Vector{Float64},
                        feas::Bool, desc::String)
    (verb == 2 || verb == 1) && println("solved $desc: ", val, " ", x, " ", feas)
    return nothing
end

for (typ, fd) in SUBPROB_SYM
    @eval function print_summary!(s::$typ, v::Int64, r::SIPSubResult, i::Int = 0)
        if i == 0
            summary_inner!(v, r.$fd.obj_val, r.$fd.sol, r.$fd.feas, "$s ")
        else
            summary_inner!(v, r.$fd.obj_val, r.$fd.sol, r.$fd.feas, "$s($i) ")
        end
    end
end

# Check convergence
function check_convergence(result::SIPResult, atol::Float64, verb::Int64)
    if abs(result.upper_bound -  result.lower_bound) < atol
        (verb == 2 || verb == 1) && println("Algorithm Converged")
        return true
    end
    return false
end

get_bnds(s::Union{LowerLevel1,LowerLevel2,LowerLevel3}, p::SIPProblem) = p.p_l, p.p_u, p.np
get_bnds(s::Union{LowerProblem,UpperProblem, ResProblem}, p::SIPProblem) = p.x_l, p.x_u, p.nx

function bnd_check(is_local::Bool, t::MOI.TerminationStatusCode, r::MOI.ResultStatusCode)
    if t === MOI.OPTIMAL
        return true
    elseif !is_local
        error("Lower problem did not solve to global optimality.
               Termination status = $t. Primal status = $r")
    elseif is_local && !((t == MOI.LOCALLY_SOLVED) || (t == MOI.ALMOST_LOCALLY_SOLVED))
        error("Lower problem did not solve to local optimality.")
    end
    return false
end
