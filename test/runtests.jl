#!/usr/bin/env julia

using Compat
using Compat.Test
using EAGO, JuMP, MathOptInterface, StaticArrays, Ipopt
const MOI = MathOptInterface

include("branch_bound.jl")
#include("ParametricInterval/runtests.jl")
include("mccormick.jl")
include("domain_reduction.jl")
include("relaxations.jl")
include("optimizer.jl")
include("script_optimizer.jl")
#include("semiinfinite.jl")
