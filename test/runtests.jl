#!/usr/bin/env julia

using Test, Printf, EAGO, MathOptInterface, Cbc, JuMP, Ipopt

const MOI = MathOptInterface
const MOIT = MOI.Test
const MOIU = MOI.Utilities
const MOIB = MOI.Bridges


#include(joinpath(@__DIR__, "moit_tests.jl"))
include(joinpath(@__DIR__, "minlp_tests.jl"))

#include(joinpath(@__DIR__, "branch_bound.jl"))
#include(joinpath(@__DIR__, "domain_reduction.jl"))
#include(joinpath(@__DIR__, "optimizer.jl"))
#include(joinpath(@__DIR__, "script_optimizer.jl"))
#include(joinpath(@__DIR__, "semiinfinite.jl"))
