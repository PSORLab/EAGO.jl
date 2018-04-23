#!/usr/bin/env julia

using EAGO

println("BEGIN TESTING INTERVAL LIBRARY...")
println("Testing Arithmetic...")
t = @elapsed include("Interval/runtests.jl")
println("done (took $t seconds).")
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

println("BEGIN TESTING NLP SOLVER")
include("NLPSolver/ExplicitOpt.jl")
println("END TESTING NLP SOLVER")

println("BEGIN TESTING SIP SOLVER")
println("END TESTING SIP SOLVER")
