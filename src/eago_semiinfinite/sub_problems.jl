# Load a model
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

# Shared LLP subroutines
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
function sipRes_llp!(t::DefaultExt, alg::A, s::S, result::SIPResult,
                     buffer::SIPBuffer, prob::SIPProblem, cb::SIPCallback,
                     i::Int64) where {A <: AbstractSIPAlgo, S <: AbstractSubproblemType}

    # build the model
    m_llp, p = build_model(t, alg, s, prob)

    # define the objective
    g(p...) = cb.gSIP[i](xbar, p)
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
function sipRes_llp!(t::ExtensionType, alg::A, s::S, result::SIPResult,
                     buffer::SIPBuffer, prob::SIPProblem, cb::SIPCallback,
                     i::Int64) where {A <: AbstractSIPAlgo, S <: AbstractSubproblemType}
    sipRes_llp!(DefaultSubproblem(), s, result, buffer, prob, cb, i)
end

# Shared bounding problems
function sipRes_bnd!(t::ExtensionType, alg::A, s::S, buffer::SIPBuffer,
                    eps_g::Float64, result::SIPResult, prob::SIPProblem,
                    cb::SIPCallback) where {A <: AbstractSIPAlgo, S <: AbstractSubproblemType}
    sipRes_bnd!(DefaultExt(), alg, s, buffer, eps_g, result, prob, cb)
end
function sipRes_bnd!(t::DefaultExt, alg::A, s::S, buffer::SIPBuffer,
                    eps_g::Float64, result::SIPResult, prob::SIPProblem,
                    cb::SIPCallback) where {A <: AbstractSIPAlgo, S <: AbstractSubproblemType}

    # create JuMP model
    m_bnd, x = build_model(t, alg, s, prob)
    disc_set = get_disc_set(s, prob)

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
    nl_obj = isone(prob.nx) ? :(($obj)(x[1])) : :(($obj)(x...))
    set_NL_objective(m_bnd, MOI.MIN_SENSE, nl_obj)

    # optimize model and check status
    JuMP.optimize!(m_bnd)
    t_status = JuMP.termination_status(m_bnd)
    r_status = JuMP.primal_status(m_bnd)
    is_feasible = bnd_check(prob.local_solver, t_status, r_status, eps_g)

    # fill buffer with subproblem result info
    buffer.is_feasible = is_feasible
    buffer.obj_value = obj_factor*JuMP.objective_bound(m_bnd)
    @__dot__ buffer.x_bar = JuMP.value(x)
    result.solution_time += MOI.get(m_bnd, MOI.SolveTime())

    return nothing
end

# get optimizer for use in subproblems
function get_sip_optimizer(t::DefaultExt, alg::A, s::S) where {A<:AbstractSIPAlgo, S<:AbstractSubproblemType}
    return EAGO.Optimizer
end
function get_sip_optimizer(t::ExtensionType, alg::A, s::AbstractSubproblemType) where {A <: AbstractSIPAlgo}
    return get_sip_optimizer(DefaultExt(), s)
end

# Printing
function print_int!(prob::SIPProblem, k::Int64, lbd::Float64, ubd::Float64,
                    eps::Float64, r::Float64, ismin::Bool)

    if (prob.verbosity == 1 || prob.verbosity == 2)

        # prints header line every hdr_intv times
        if (mod(k, prob.hdr_intv) == 0 || k == 1)
            println("| Iteration | Lower Bound | Upper Bound |   eps   |   r   |  Gap  |  Ratio  |")
        end

        # prints iteration summary every prnt_intv times
        if mod(k, prob.prt_intv) == 0

            print_str = "| "

            max_len = 15
            temp_str = string(k)
            len_str = length(temp_str)
            print_str *= (" "^(max_len - len_str))*temp_str*" | "

            max_len = 15
            lower_adj = ismin ? lbd : ubd
            temp_str = string(round(lower_adj, digits = 6))
            len_str = length(temp_str)
            print_str *= (" "^(max_len - len_str))*temp_str*" | "

            max_len = 15
            upper_adj = ismin ? ubd : lbd
            temp_str = string(round(upper_adj, digits = 6))
            len_str = length(temp_str)
            print_str *= (" "^(max_len - len_str))*temp_str*" | "

            max_len = 15
            temp_str = string(round(eps, digits = 6))
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

function print_summary!(verb::Int64, val::Float64, x::Vector{Float64},
                          feas::Bool, desc::String)
    (verb == 2 || verb == 1) && println("solved $desc: ", val, " ", x, " ", feas)
    return nothing
end

# Check convergence
function check_convergence(LBD::Float64, UBD::Float64, atol::Float64, verb::Int64)
    if abs(UBD - LBD) < atol
        (verb == 2 || verb == 1) && println("Algorithm Converged")
        return true
    end
    return false
end


get_bnds(s::Union{LowerLevel1,LowerLevel2}, p::SIPProblem) = p.pL, p.pU. p.np
get_bnds(s::Union{LowerProblem,UpperProblem}, p::SIPProblem) = p.xL, p.xU. p.nx

get_sip_optimizer(t::DefaultExt, alg::SIPRes, s::S) where S <: AbstractSubproblemType = EAGO.Optimizer

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
