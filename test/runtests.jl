#!/usr/bin/env julia

using Test
using EAGO, JuMP, MathOptInterface, StaticArrays, Ipopt
using EAGO.McCormick
using IntervalArithmetic
const MOI = MathOptInterface

include("branch_bound.jl")
include("mccormick.jl")
include("domain_reduction.jl")
#include("relaxations.jl")
include("optimizer.jl")
include("script_optimizer.jl")
#include("semiinfinite.jl")
