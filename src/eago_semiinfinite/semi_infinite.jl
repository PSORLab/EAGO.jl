# Copyright (c) 2018: Matthew Wilhelm & Matthew Stuber.
# This code is licensed under MIT license (see LICENSE.md for full details)
#############################################################################
# EAGO
# A development environment for robust and global optimization
# See https://github.com/PSORLab/EAGO.jl
#############################################################################
# src/eago_semiinfinite/semi_infinite.jl
# Defines high level interface for SIP routines.
#############################################################################

include("types.jl")
include("subproblems.jl")

function sip_solve!(t::AbstractSIPAlgo, subresult, prob, result, cb)
    error("Algorithm $t not supported by explicit_sip_solve.")
end

"""
    sip_solve

Solve an SIP with decision variable bounds `x_l` to `x_u`, uncertain variable
bounds `p_l` to `p_u`, an objective function of `f`, and `gSIP` seminfiniite
constraint(s).
"""
function sip_solve(t::T, x_l::Vector{Float64}, x_u::Vector{Float64},
                   p_l::Vector{Float64}, p_u::Vector{Float64},
                   f::Function, gSIP; kwargs...) where T <: AbstractSIPAlgo

    @assert length(p_l) == length(p_u)
    @assert length(x_l) == length(x_u)

    prob = SIPProblem(x_l, x_u, p_l, p_u, gSIP, m, kwargs)
    subresult = SIPSubResult(prob.np, prob.nx)
    result = SIPResult(prob.nx, prob.np)
    cb = SIPCallback(f, gSIP)

    sip_solve!(t, subresult, prob, result, cb)
    return result
end

function sip_solve(t::T, x_l::Vector{Float64}, x_u::Vector{Float64},
                   p_l::Vector{Float64}, p_u::Vector{Float64},
                   f::Function, gSIP::Function; kwargs...) where T <: AbstractSIPAlgo
    return sip_solve(t, x_l, x_u, p_l, p_u, f, [gSIP], kwargs...)
end

include("algorithms/sip_hybrid.jl")
include("algorithms/sip_res_rev.jl")
include("algorithms/sip_res.jl")
