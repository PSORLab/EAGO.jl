# Copyright (c) 2018 Matthew Wilhelm, Robert Gottlieb, Dimitri Alston,
# Alireza Miraliakbar, Matthew Stuber, and the University of Connecticut (UConn)
# This code is licensed under the MIT license (see LICENSE.md for full details).
################################################################################
# EAGO
# A development environment for robust and global optimization.
# https://github.com/PSORLab/EAGO.jl
################################################################################
# src/eago_semiinfinite/semi_infinite.jl
# Defines a high-level interface for SIP routines.
################################################################################

include("types.jl")
include("subproblems.jl")

function sip_solve(t::ExtensionType, alg::AbstractSIPAlgo, subresult, prob, result, cb)
    error("Algorithm $t not supported by explicit_sip_solve.")
end

"""
    sip_solve

Solve an SIP with decision variable bounds `x_l` to `x_u`, uncertain variable
bounds `p_l` to `p_u`, an objective function of `f`, and `gSIP` seminfiniite
constraint(s).
"""
function sip_solve(alg::T, x_l::Vector{Float64}, x_u::Vector{Float64},
                   p_l::Vector{Float64}, p_u::Vector{Float64},
                   f::Function, gSIP::Vector; d = DefaultExt(), kwargs...) where {T <: AbstractSIPAlgo}

    @assert length(p_l) == length(p_u)
    @assert length(x_l) == length(x_u)

    prob = SIPProblem(x_l, x_u, p_l, p_u, gSIP, kwargs)
    subresult = SIPSubResult(prob.nx, prob.np, prob.nSIP, prob.abs_tolerance)
    result = SIPResult(prob.nx, prob.np)
    cb = SIPCallback(f, gSIP)

    sip_solve!(d, alg, subresult, prob, result, cb)
    return result
end

function sip_solve(alg::T, x_l::Vector{Float64}, x_u::Vector{Float64},
                   p_l::Vector{Float64}, p_u::Vector{Float64},
                   f::Function, gSIP::Function; d = DefaultExt(), kwargs...) where {T <: AbstractSIPAlgo}
    return sip_solve(alg, x_l, x_u, p_l, p_u, f, [gSIP], d = d, kwargs...)
end

include("nonconvex_algorithms/sip_hybrid.jl")
include("nonconvex_algorithms/sip_res_rev.jl")
include("nonconvex_algorithms/sip_res.jl")
