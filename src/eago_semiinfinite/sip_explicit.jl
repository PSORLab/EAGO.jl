# Copyright (c) 2018: Matthew Wilhelm & Matthew Stuber.
# This code is licensed under MIT license (see LICENSE.md for full details)
#############################################################################
# EAGO
# A development environment for robust and global optimization
# See https://github.com/PSORLab/EAGO.jl
#############################################################################
# src/eago_semiinfinite/sip_explicit.jl
# An implementation of the SIPres algorithm.
#############################################################################

abstract type AbstractSubproblemType end
struct LowerLevel1 <: AbstractSubproblemType end
struct LowerLevel2 <: AbstractSubproblemType end
struct LowerProblem <: AbstractSubproblemType end
struct UpperProblem <: AbstractSubproblemType end

get_sip_optimizer(t::DefaultExt, alg::SIPRes, s::S) where S <: AbstractSubproblemType = EAGO.Optimizer
get_sip_optimizer(t::ExtensionType, alg::SIPRes, s::AbstractSubproblemType) = get_sip_optimizer(DefaultExt(), s)

get_sip_kwargs(alg::SIPRes, s::LowerLevel1, p::SIPProblem) = p.kwargs_llp1
get_sip_kwargs(alg::SIPRes, s::LowerLevel2, p::SIPProblem) = p.kwargs_llp2
get_sip_kwargs(alg::SIPRes, s::LowerProblem, p::SIPProblem) = p.kwargs_lbd
get_sip_kwargs(alg::SIPRes, s::UpperProblem, p::SIPProblem) = p.kwargs_ubd

get_bnds(s::Union{LowerLevel1,LowerLevel2}, p::SIPProblem) = p.pL, p.pU. p.np
get_bnds(s::Union{LowerProblem,UpperProblem}, p::SIPProblem) = p.xL, p.xU. p.nx

# Load a model in a way that gets rid of any issues with pointers for C references
# particularly in the EAGO solver...
function build_model(t::DefaultExt, a::A, s::S, p::SIPProblem) where {A <: AbstractSIPAlgo, S <: AbstractSubproblemType}
    model = Model(get_sip_optimizer(t,a,s))
    for (k,v) in get_sip_kwargs(a,s,p)
        MOI.set(model, MOI.RawParameter(String(k)), v)
    end
    vL, vU, nv = get_bnds(s,a,p)
    @variable(model, vL[i] <= v[i=1:nv] <= vU[i])
    return model, v
end
function build_model(t::ExtensionType, a::A, s::S, p::SIPProblem) where {A <: AbstractSIPAlgo, S <: AbstractSubproblemType}
    build_model(DefaultExt(), s, p)
end

function add_uncertainty_constraint!(model::JuMP.Model, problem::SIPProblem)
    #if !isnothing(model.polyhedral_uncertainty_set)
        #@constraint(model, )
    #end
    #if !isnothing(model.ellipsodial_uncertainty_set)
        #@constraint(model, )
    #end
    return nothing
end

function llp_check(islocal::Bool, t::MOI.TerminationStatusCode, r::MOI.PrimalResultCode)
    valid, feasible = is_globally_optimal(t, r)
    if islocal && ((t != MOI.LOCALLY_SOLVED) && (t != MOI.ALMOST_LOCALLY_SOLVED))
        error("Lower problem did not solve to local optimality.")
    elseif !valid_result
        error("Error in lower level problem. Termination status = $t, primal status = $r.")
    end
    return is_feasible
end

function sipRes_llp!(t::DefaultExt, alg::SIPRes, s::S, result::SIPResult,
                     buffer::SIPBuffer, prob::SIPProblem, cb::SIPCallback,
                     indx::Int64) where S <: AbstractSubproblemType

    # build the model
    m_llp, p = build_model_llp(t, alg, s, prob)

    # define the objective
    g(p...) = cb.gSIP[indx](xbar, p)
    register(m_llp, :g, prob.np, g, autodiff=true)
    nl_obj = isone(prob.np) ? :(-($g)(p[1])) :  :(-($g)(p...))
    set_NL_objective(m_llp, MOI.MIN_SENSE, nl_obj)

    # add uncertainty constraints
    add_uncertainty_constraint!(m_llp, problem)

    # optimize model and check status
    JuMP.optimize!(m_llp)
    tstatus = JuMP.termination_status(m_llp)
    rstatus = JuMP.primal_status(m_llp)
    is_feasible = llp_check(prob.local_solver, tstatus, rstatus)

    # fill buffer with subproblem result info
    buffer.is_feasible = is_feasible
    buffer.obj_value = -JuMP.objective_value(m_llp)
    @__dot__ buffer.p_bar = JuMP.value(p)
    result.solution_time += MOI.get(m_llp, MOI.SolveTime())

    return nothing
end
function sipRes_llp!(t::ExtensionType, alg::SIPRes, s::S, result::SIPResult,
                     buffer::SIPBuffer, prob::SIPProblem, cb::SIPCallback,
                     indx::Int64) where S <: AbstractSubproblemType
    sipRes_llp!(DefaultSubproblem(), als, s, result, buffer, prob, cb, indx)
end

function bnd_check(is_local::Bool, t::MOI.TerminationStatusCode,
                   r::MOI.PrimalResultCode, eps_g::Float64)
    valid_result, is_feasible = is_globally_optimal(t, r)
    if (!(valid_result && is_feasible) && iszero(eps_g)) && !is_local
        error("Lower problem did not solve to global optimality.
               Termination status = $t. Primal status = $r")
    elseif (!(valid_result && is_feasible) && iszero(eps_g)) && is_local &&
           !((t == MOI.LOCALLY_SOLVED) || (t == MOI.ALMOST_LOCALLY_SOLVED))
        error("Lower problem did not solve to local optimality.")
    end
    return is_feasible
end

function sipRes_bnd(t::AbstractSubproblemType, buffer, initialize_extras, disc_set::Vector{Vector{Vector{Float64}}},
                    eps_g::Float64, sip_storage::SIPResult, prob::SIPProblem, flag::Bool, cb::SIPCallback)
    sipRes_bnd(t, buffer, initialize_extras, disc_set, eps_g, sip_storage, prob, flag, cb)
end

# should be done
function sipRes_bnd(t::DefaultSubproblem, buffer, initialize_extras, disc_set::Vector{Vector{Vector{Float64}}},
                    eps_g::Float64, sip_storage::SIPResult,
                    prob::SIPProblem, flag::Bool, cb::SIPCallback)

    # create JuMP model
    m_bnd = build_model(prob)
    @variable(m_bnd, prob.x_l[i] <= x[i=1:prob.nx] <= prob.x_u[i])
    initialize_extras(m_bnd, x)                    # TODO: Still worth having this feature?

    for i = 1:prob.nSIP
        for j = 1:length(disc_set)
            gi = Symbol("g$i$j")
            g(x...) = cb.gSIP[i](x, disc_set[j][i])
            register(m_bnd, gi, nx, g, autodiff=true)
            xids = VariableRef[VariableRef(m_bnd, MOI.VI(i)) for i = 1:prob.nx]
            JuMP.add_NL_constraint(m_bnd, :($(gi)($(tuple(xids...))) + $ep_g <= 0))
        end
    end

    # define the objective
    obj_factor = prob.sense == :min ? 1.0 : -1.0
    obj(x...) =  obj_factor*cb.f(x)
    register(m_bnd, :obj, prob.nx, obj, autodiff=true)
    nl_obj = isone(prob.nx) ? :(($obj)(x[1])) :(($obj)(x...))
    set_NL_objective(m_bnd, MOI.MIN_SENSE, nl_obj)

    # optimize model and check status
    JuMP.optimize!(m_bnd)
    termination_status = JuMP.termination_status(m_bnd)
    result_status = JuMP.primal_status(m_bnd)
    is_feasible = bnd_check(prob.local_solver, termination_status, result_status, eps_g)

    # fill buffer with subproblem result info
    buffer.is_feasible = is_feasible
    buffer.obj_value = obj_factor*JuMP.objective_bound(m_bnd)
    @__dot__ buffer.x_bar = JuMP.value(x)
    result.solution_time += MOI.get(model_llp, MOI.SolveTime())

    return false
end

function sip_solve!(t::AbstractSIPAlgo, init_bnd, prob, result, cb)
    error("Algorithm $t not supported by explicit_sip_solve.")
end

function sip_solve!(t::SIPRes, init_bnd, prob::SIPProblem, result::SIPResult, cb::SIPCallback)

    verb = prob.verbosity
    header_interval = prob.header_interval
    print_interval = prob.print_interval
    abs_tolerance = prob.absolute_tolerance

    # initializes solution
    xbar = (prob.x_u + prob.x_l)/2.0
    pbar = (prob.p_u + prob.p_l)/2.0
    nx = prob.nx
    np = prob.np
    tolerance = prob.constraint_tolerance
    ismin = prob.sense === :min

    eps_g =  prob.initial_eps_g
    r = prob.initial_r
    feas = true
    xstar = fill(NaN,(nx,))
    new_disc_points = fill(zeros(Float64, np), )

    # checks for convergence
    for k = 1:prob.iteration_limit

        # check for termination
        check_convergence(result.lower_bound, result.upper_bound, abs_tolerance, verb) && (break)

        # solve lower bounding problem and check feasibility
        sipRes_bnd!(buffer, init_bnd, prob.lower_disc, 0.0, result, prob, true, cb)
        result.lower_bound = buffer.objective_value
        if !is_feasible(buffer)
            result.feasibility = feas
            println("Lower Bounding Problem Not Feasible. Algorithm Terminated")
            break
        end
        print_summary!(verb, result.lower_bound, result.xsol, result.feasibility, "LBD")

        # solve inner program  and update lower discretization set
        non_positive_flag = true
        for i = 1:prob.nSIP
            sipRes_llp!(t, alg, LowerLevel1(), result, buffer, prob, cb, i)
            if (verb == 1) || (verb == 2)
                print_summary!(buffer, "Lower LLP$i")
            end
            if (llp_out[1] + tolerance > 0.0)
                non_positive_flag = false
            end
            buffer.lower_disc[i] .= buffer.pbar
        end
        push!(prob.lower_disc, deepcopy(buffer.lower_disc))
        if non_positive_flag
            result.feasibility = true
            break
        end

        # solve upper bounding problem, if feasible solve lower level problem,
        # and potentially update upper discretization set
        sipRes_bnd!(buffer, init_bnd, prob.upper_disc, eps_g, result, prob, false, cb)
        print_summary!(verb, buffer, :x, "UBD")
        if is_feasible(buffer)
            non_positive_flag = true
            for i = 1:prob.nSIP
                sipRes_llp!(t, alg, LowerLevel2(), result, buffer, prob, cb, i)
                print_summary!(verb, buffer, :p, "Upper LLP$i")
                if llp2_out[1] + tolerance/10.0 > 0.0
                    non_positive_flag = false
                end
                buffer.upper_disc[i] .= buffer.pbar
            end
            if non_positive_flag
                if (val <= result.upper_bound) && ismin
                    result.upper_bound = val
                    xstar .= result.xsol
                elseif (val <= result.upper_bound) && ~ismin && (result.upper_bound === Inf)
                    result.upper_bound = val
                    xstar .= result.xsol
                elseif (val >= result.upper_bound) && ~ismin
                    result.upper_bound = val
                    xstar .= result.xsol
                end
                eps_g /= r
            else
                push!(prob.upper_disc, deepcopy(buffer.upper_disc))
            end
        else
            eps_g /= r
        end

        # print iteration information and advance
        print_int!(verb, header_interval, print_interval, k, result.lower_bound, result.upper_bound, eps_g, r, ismin)
        result.iteration_number = k
    end

    return nothing
end

"""
  explicit_sip_solve

Solve an SIP with decision variable bounds `x_l` to `x_u`, uncertain variable
bounds `p_l` to `p_u`, an objective function of `f`, and `gSIP` seminfiniite
constraint(s).
"""
function sip_solve(t::T, x_l::Vector{Float64}, x_u::Vector{Float64},
                   p_l::Vector{Float64}, p_u::Vector{Float64},
                   f::Function, gSIP; kwargs...) where T <: AbstractSIPAlgo

    @assert length(p_l) == length(p_u)
    @assert length(x_l) == length(x_u)

    # collects all keyword arguments of the form :gSIP1, gSIP24234, :gSIPHello
    m = haskey(kwargs, :sip_optimizer) ? kwargs[:sip_optimizer] : EAGO.Optimizer       # TODO: Optimizer type is runtime variable... need to key to type...

    prob = SIPProblem(x_l, x_u, p_l, p_u, gSIP, m, kwargs)
    buffer = SIPBuffer(prob.np, prob.nx)
    result = SIPResult(prob.nx, prob.np)
    cb = SIPCallback(f, gSIP)

    sip_solve!(t, init_bnd, prob, result, cb)
    return result
end

function sip_solve(t::T, x_l::Vector{Float64}, x_u::Vector{Float64},
                            p_l::Vector{Float64}, p_u::Vector{Float64},
                            f::Function, gSIP::Function; kwargs...) where T <: AbstractSIPAlgo
    return sip_solve(t, x_l, x_u, p_l, p_u, f, [gSIP], kwargs...)
end
