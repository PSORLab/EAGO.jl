#!/usr/bin/env julia

using Test, Printf, EAGO, MathOptInterface, Cbc, JuMP

const MOI = MathOptInterface
const MOIT = MOI.Test
const MOIU = MOI.Utilities
const MOIB = MOI.Bridges


include(joinpath(@__DIR__, "moit_tests.jl"))
include(joinpath(@__DIR__, "minlp_tests.jl"))

#include("branch_bound.jl")
#include("domain_reduction.jl")
#include("optimizer.jl")
#include("script_optimizer.jl")
#include("semiinfinite.jl")
