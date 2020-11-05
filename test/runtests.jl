#!/usr/bin/env julia

using Test
using EAGO, JuMP, MathOptInterface, Ipopt, ForwardDiff
using IntervalArithmetic, SpecialFunctions
const MOI = MathOptInterface
using ForwardDiff: Dual, Partials

include("branch_bound.jl")
include("domain_reduction.jl")
include("optimizer.jl")
include("script_optimizer.jl")
include("semiinfinite.jl")
