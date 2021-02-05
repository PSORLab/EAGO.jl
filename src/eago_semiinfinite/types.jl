# Copyright (c) 2018: Matthew Wilhelm & Matthew Stuber.
# This code is licensed under MIT license (see LICENSE.md for full details)
#############################################################################
# EAGO
# A development environment for robust and global optimization
# See https://github.com/PSORLab/EAGO.jl
#############################################################################
# # src/eago_semiinfinite/types.jl
# Defines intermediate types used by SIP routines.
#############################################################################

abstract type AbstractSIPAlgo end

abstract type AbstractSubproblemType end
struct LowerLevel1 <: AbstractSubproblemType end
struct LowerLevel2 <: AbstractSubproblemType end
struct LowerLevel3 <: AbstractSubproblemType end
struct LowerProblem <: AbstractSubproblemType end
struct UpperProblem <: AbstractSubproblemType end
struct ResProblem <: AbstractSubproblemType end

"""
    SIPResult

Structure storing the results of the SIPres algorithm.
"""
mutable struct SIPResult
    iteration_number::Int64
    res_iteration_number::Int64
    upper_bound::Float64
    lower_bound::Float64
    feasibility::Bool
    xsol::Vector{Float64}
    psol::Vector{Float64}
    solution_time::Float64
end
SIPResult() = SIPResult(1, 1, Inf, -Inf, true, Float64[], Float64[], 0.0)
SIPResult(nx::Int, np::Int) = SIPResult(1, Inf, -Inf, true, zeros(nx), zeros(np), 0.0)


"""
      SIPProblem

Structure storing problem information for the solution routine.
"""
mutable struct SIPProblem
    x_l::Vector{Float64}
    x_u::Vector{Float64}
    p_l::Vector{Float64}
    p_u::Vector{Float64}

    np::Int
    nSIP::Int
    nx::Int

    absolute_tolerance::Float64
    constraint_tolerance::Float64
    iteration_limit::Int
    res_iteration_limit::Int
    initial_eps_g::Float64
    initial_r::Float64

    return_hist::Bool
    header_interval::Int
    print_interval::Int
    verbosity::Int

    local_solver::Bool

    #polyhedral_uncertainty_set
    #ellipsodial_uncertainty_set

    optimizer
    kwargs
end

get_sip_kwargs(s::LowerLevel1, p::SIPProblem) = p.kwargs_llp1
get_sip_kwargs(s::LowerLevel2, p::SIPProblem) = p.kwargs_llp2
get_sip_kwargs(s::LowerProblem, p::SIPProblem) = p.kwargs_lbd
get_sip_kwargs(s::UpperProblem, p::SIPProblem) = p.kwargs_ubd
get_sip_kwargs(s::ResProblem, p::SIPProblem) = p.kwargs_res

function SIPProblem(x_l::Vector{Float64}, x_u::Vector{Float64},
                    p_l::Vector{Float64}, p_u::Vector{Float64},
                    gSIP, optimizer, kwargs)

    initial_eps_g = haskey(kwargs, :sip_initial_eps_g) ? kwargs[:sip_initial_eps_g] : 1.0
    initial_r = haskey(kwargs, :sip_initial_r) ? kwargs[:sip_initial_r] : 2.0

    (initial_r <= 1.0) && error("initial_r must be greater than 1")
    (initial_eps_g <= 0.0) && error("eps_g must be greater than 0")

    absolute_tolerance = haskey(kwargs, :sip_absolute_tolerance) ? kwargs[:sip_absolute_tolerance] : 1E-3
    constraint_tolerance = haskey(kwargs, :sip_constraint_tolerance) ? kwargs[:sip_constraint_tolerance] : 1E-3
    iteration_limit = haskey(kwargs, :sip_iteration_limit) ? kwargs[:sip_iteration_limit] : 100
    res_iteration_limit = haskey(kwargs, :sip_res_iteration_limit) ? kwargs[:sip_res_iteration_limit] : 100
    return_hist = haskey(kwargs, :sip_return_hist) ? kwargs[:sip_return_hist] : false
    header_interval = haskey(kwargs, :sip_header_interval) ? kwargs[:sip_header_interval] : 20
    print_interval = haskey(kwargs, :sip_print_interval) ? kwargs[:sip_print_interval] : 1
    verbosity = haskey(kwargs, :sip_verbosity) ? kwargs[:sip_verbosity] : 1
    local_solver = haskey(kwargs, :sip_local_solver) ?  kwargs[:sip_local_solver] : false

    np = length(p_l)
    nx = length(x_l)

    opt_dict = Dict{Symbol,Any}()
    for key in keys(kwargs)
        string_key = String(key)
        if string_key[1:3] !== "sip"
            opt_dict[key] = kwargs[key]
        end
    end

    # polyhedral_uncertainty_set = nothing
    # ellipsodial_uncertainty_set = nothing
    # conic_uncertainty_set = nothing
    # convex_uncertainty_set = nothing

    nSIP = length(gSIP)

    SIPProblem(x_l, x_u, p_l, p_u, np, nSIP, nx,
               absolute_tolerance, constraint_tolerance,
               iteration_limit,
               initial_eps_g, initial_r, return_hist, header_interval,
               print_interval, verbosity, local_solver,
               #polyhedral_uncertainty_set,
               #ellipsodial_uncertainty_set, conic_uncertainty_set,
               #convex_uncertainty_set,
               optimizer, opt_dict)
end

struct SIPCallback
    f
    gSIP
end

mutable struct SubProblemInfo
    sol::Vector{Float64}
    obj_val::Float64
    obj_bnd::Float64
    feas::Bool
    tol::Vector{Float64}
end
function SubProblemInfo(nd::Int, ng::Int, tol::Float64)
    SubProblemInfo(zeros(nd), 0.0, 0.0, false, fill(tol, ng))
end

"""
    SIPBuffer

Hold objective value, solution, discretization set, and feasibility status of
each subproblem encountered by SIP algorithm.
"""
Base.@kwdef mutable struct SIPSubResult
    fRes::Float64 = Inf
    lbd::SubProblemInfo
    ubd::SubProblemInfo
    res::SubProblemInfo
    llp1::SubProblemInfo
    llp2::SubProblemInfo
    llp3::SubProblemInfo
    r_g::Float64 = 5.0
    r_l::Float64 = 5.0
    r_u::Float64 = 5.0
    eps_g::Vector{Float64} = Float64[]
    eps_l::Vector{Float64} = Float64[]
    eps_u::Vector{Float64} = Float64[]
    disc_l_buffer::Vector{Float64} = Float64[]
    disc_u_buffer::Vector{Float64} = Float64[]
    disc_l::Vector{Vector{Vector{Float64}}} = Vector{Vector{Float64}}[]
    disc_u::Vector{Vector{Vector{Float64}}} = Vector{Vector{Float64}}[]
end
function SIPSubResult(nx::Int, np::Int, ng::Int, tol::Float64)
    buffer = SIPSubResult(lbd  = SubProblemInfo(nx, 1, tol),
                          ubd  = SubProblemInfo(nx, 1, tol),
                          res  = SubProblemInfo(nx, 1, tol),
                          llp1 = SubProblemInfo(np, ng, tol),
                          llp2 = SubProblemInfo(np, ng, tol),
                          llp3 = SubProblemInfo(np, ng, tol))
    append!(buffer.eps_g, fill(1E-3, ng))
    append!(buffer.eps_l, fill(1E-3, ng))
    append!(buffer.eps_u, fill(1E-3, ng))
    append!(buffer.disc_l_buffer, zeros(np))
    append!(buffer.disc_u_buffer, zeros(np))
    return buffer
end

get_eps(::LowerProblem, sr::SIPSubResult, i::Int) = 0.0
get_eps(::UpperProblem, sr::SIPSubResult, i::Int) = sr.eps_g[i]

const SUBPROB_SYM = Dict{Symbol,Symbol}(:LowerProblem => :lbd,
                                        :UpperProblem => :ubd,
                                        :ResProblem   => :res,
                                        :LowerLevel1  => :llp1,
                                        :LowerLevel2  => :llp1,
                                        :LowerLevel3  => :llp3)

for (typ, fd) in SUBPROB_SYM
    @eval function load!(::$typ, sr::SIPSubResult, feas::Bool, val::Float64,
                                 bnd::Float64, x::Vector{Float64})
        sr.$fd.feas = feas
        sr.$fd.obj_val = obj
        sr.$fd.obj_bnd = bnd
        sr.$fd.sol .= x
        return nothing
    end
end
