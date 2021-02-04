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

function set_tolerance!(t::DefaultExt, alg::AbstractSIPAlgo, s::S, m::JuMP.Model,
                        sr::SIPSubResult, i::Int) where {S <: AbstractSubproblemType}
    return nothing
end
function set_tolerance!(t::ExtensionType, alg::A, s::S, m::JuMP.Model,
                        sr::SIPSubResult, i::Int) where {A <: AbstractSIPAlgo,
                                                         S <: AbstractSubproblemType}
    set_tolerance!(t, alg, s, m, sr, i)
end

function get_disc_set(t::DefaultExt, alg::AbstractSIPAlgo, s::S, sr::SIPProblem, i::Int) where {S <: AbstractSubproblemType}
    return Vector{Float64}[]
end
function get_disc_set(t::ExtensionType, alg::AbstractSIPAlgo, s::S, sr::SIPProblem, i::Int) where {S <: AbstractSubproblemType}
    get_disc_set(DefaultExt(), alg, s, p, i)
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
function sip_llp!(t::DefaultExt, alg::A, s::S, result::SIPResult,
                  sr::SIPSubResult, prob::SIPProblem, cb::SIPCallback,
                  i::Int64, tol::Float64 = -Inf) where {A <: AbstractSIPAlgo, S <: AbstractSubproblemType}

    # build the model
    m, p = build_model(t, alg, s, prob)
    set_tolerance!(t, alg, s, m, sr, i)

    # define the objective
    g(p...) = cb.gSIP[i](xbar, p)
    register(m, :g, prob.np, g, autodiff=true)
    nl_obj = isone(prob.np) ? :(-($g)(p[1])) :  :(-($g)(p...))
    set_NL_objective(m, MOI.MIN_SENSE, nl_obj)

    # add uncertainty constraints
    add_uncertainty_constraint!(m, problem)

    # optimize model and check status
    JuMP.optimize!(m)
    tstatus = JuMP.termination_status(m)
    rstatus = JuMP.primal_status(m)
    feas = llp_check(prob.local_solver, tstatus, rstatus)

    # fill buffer with subproblem result info
    load!(s, buffer, feas, -JuMP.objective_bound(m), JuMP.value(x))
    result.solution_time += MOI.get(m, MOI.SolveTime())

    return nothing
end
function sip_llp!(t::ExtensionType, alg::A, s::S, result::SIPResult,
                     sr::SIPSubResult, prob::SIPProblem, cb::SIPCallback,
                     i::Int64) where {A <: AbstractSIPAlgo, S <: AbstractSubproblemType}
    sip_llp!(DefaultSubproblem(), s, result, buffer, prob, cb, i)
end

# Shared bounding problems
function sip_bnd!(t::ExtensionType, alg::A, s::S, sr::SIPSubResult, result::SIPResult,
                  prob::SIPProblem, cb::SIPCallback) where {A <: AbstractSIPAlgo,
                                                            S <: AbstractSubproblemType}
    sip_bnd!(DefaultExt(), alg, s, buffer, eps_g, result, prob, cb)
end
function sip_bnd!(t::DefaultExt, alg::A, s::S, sr::SIPSubResult, result::SIPResult,
                  prob::SIPProblem, cb::SIPCallback) where {A <: AbstractSIPAlgo,
                                                            S <: AbstractSubproblemType}

    # create JuMP model
    m, x = build_model(t, alg, s, prob)

    for i = 1:prob.nSIP
        ϵ_g = get_eps(s, sr, i)
        disc_set = get_disc_set(t, alg, s, sr, i)
        for j = 1:length(disc_set)
            gi = Symbol("g$i$j")
            g(x...) = cb.gSIP[i](x, disc_set[i][j])
            register(m, gi, nx, g, autodiff=true)
            JuMP.add_NL_constraint(m, :($(gi)($(tuple(x...))) + $ϵ_g <= 0))
        end
    end

    # define the objective
    obj(x...) = cb.f(x)
    register(m, :obj, prob.nx, obj, autodiff=true)
    nl_obj = isone(prob.nx) ? :(($obj)(x[1])) : :(($obj)(x...))
    set_NL_objective(m, MOI.MIN_SENSE, nl_obj)

    # optimize model and check status
    JuMP.optimize!(m)
    t_status = JuMP.termination_status(m)
    r_status = JuMP.primal_status(m)
    feas = bnd_check(prob.local_solver, t_status, r_status, eps_g)

    # fill buffer with subproblem result info
    load!(s, buffer, feas, JuMP.objective_bound(m), JuMP.value(x))
    result.solution_time += MOI.get(m, MOI.SolveTime())

    return nothing
end

# Adaptive restriction subproblem
function sip_res!(t::ExtensionType, alg::A, s::S, sr::SIPSubResult,
                    eps_g::Float64, result::SIPResult, prob::SIPProblem,
                    cb::SIPCallback) where {A <: AbstractSIPAlgo, S <: AbstractSubproblemType}
    sip_res!(DefaultExt(), alg, s, buffer, eps_g, result, prob, cb)
end
function sip_res!(t::DefaultExt, alg::A, s::S, sr::SIPSubResult,
                  eps_g::Float64, result::SIPResult, prob::SIPProblem,
                  cb::SIPCallback) where {A <: AbstractSIPAlgo, S <: AbstractSubproblemType}

    # create JuMP model & variables
    m, x = build_model(t, alg, s, prob)
    @variable(m, η)

    # add discretized semi-infinite constraint
    for i = 1:prob.nSIP
        disc_set = get_disc_set(t, alg, s, sr, i)
        for j = 1:length(disc_set)
            gi = Symbol("g$i$j")
            g(x...) = cb.gSIP[i](x, disc_set[i][j])
            register(m, gi, nx, g, autodiff=true)
            JuMP.add_NL_constraint(m, :($(gi)($(tuple(x...))) + η <= 0))
        end
    end

    # add epigraph reformulated objective
    register(m, :f, nx, cb.f, autodiff=true)
    fRes = #TODO
    JuMP.add_NL_constraint(m, :(f($(tuple(x...))) + $fRes <= 0))

    # define the objective
    @objective(m, Min, -η)

    # optimize model and check status
    JuMP.optimize!(m)
    t_status = JuMP.termination_status(m)
    r_status = JuMP.primal_status(m)
    feas = llp_check(prob.local_solver, t_status, r_status)

    # fill buffer with subproblem result info
    load!(s, buffer, feas, JuMP.objective_bound(m), JuMP.value(x))
    result.solution_time += MOI.get(m, MOI.SolveTime())

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
function print_int!(verb::Int, prob::SIPProblem, k::Int64, result::SIPResult,
                    eps::Float64, r::Float64)

    lbd = result.lower_bound
    ubd = result.upper_bound

    if (prob.verbosity == 1 || prob.verbosity == 2)

        # prints header line every hdr_intv times
        if (mod(k, prob.hdr_intv) == 0 || k == 1)
            println("| Iteration | Lower Bound | Upper Bound |   r   |  Gap  |  Ratio  |")
        end

        # prints iteration summary every prnt_intv times
        if mod(k, prob.prt_intv) == 0

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
    function print_summary!(s::$typ, v::Int64, r::SIPSubResult, i::Int = 0)
        if i == 0
            summary_inner!(v, r.$fd.obj_value, r.$fd.sol, r.$fd.feas, "$s ")
        else
            summary_inner!(v, r.$fd.obj_value, r.$fd.sol, r.$fd.feas, "$s($i) ")
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
