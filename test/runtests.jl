#!/usr/bin/env julia

using EAGO

println("BEGIN TESTING INTERVAL LIBRARY...")
include("Interval/runtests.jl")
println("TEST INTERVAL LIBRARY COMPLETE.")

println("BEGIN TESTING BRANCH AND BOUND LIBRARY...")
include("BranchBound/runtests.jl")
println("TESTING BRANCH AND BOUND LIBRARY COMPLETE.")

println("BEGIN TESTING PARAMETRIC INTERVAL LIBRARY...")
include("ParametricInterval/runtests.jl")
println("TESTING PARAMETRIC INTERVAL LIBRARY COMPLETE.")

println("BEGIN TESTING MCCORMICK LIBRARY...")
include("McCormick/runtests.jl")
println("TESTING MCCORMICK LIBRARY COMPLETE.")

println("BEGIN DOMAIN REDUCTION LIBRARY...")
include("DomainReduction/runtests.jl")
println("TESTING DOMAIN REDUCTION COMPLETE.")
