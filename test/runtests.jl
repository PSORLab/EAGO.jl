#!/usr/bin/env julia

using Compat
using Compat.Test
using EAGO, JuMP, MathOptInterface, StaticArrays, Ipopt
const MOI = MathOptInterface

println("BEGIN TESTING BRANCH AND BOUND LIBRARY...")
include("branch_bound.jl")
println("TESTING BRANCH AND BOUND LIBRARY COMPLETE.")

println("BEGIN TESTING PARAMETRIC INTERVAL LIBRARY...")
#include("ParametricInterval/runtests.jl")
println("TESTING PARAMETRIC INTERVAL LIBRARY COMPLETE.")

println("BEGIN TESTING MCCORMICK LIBRARY...")
#include("mccormick.jl")
println("TESTING MCCORMICK LIBRARY COMPLETE.")

println("BEGIN DOMAIN REDUCTION LIBRARY...")
include("domain_reduction.jl")
println("TESTING DOMAIN REDUCTION COMPLETE.")

println("BEGIN RELAXATION ROUTINES...")
include("relaxations.jl")
println("TESTING RELAXATION ROUTINES COMPLETE.")

println("BEGIN TESTING OPTIMIZER")
include("Optimizer/optimizer.jl")
println("END TESTING OPTIMIZER")

println("BEGIN TESTING SCRIPT BRIDGE")
include("script_optimizer.jl")
println("END TESTING SCRIPT BRIDGE")

println("BEGIN TESTING SIP SOLVER")
#include("semiinfinite.jl")
println("END TESTING SIP SOLVER")
