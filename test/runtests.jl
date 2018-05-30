#!/usr/bin/env julia

workspace()

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

println("BEGIN TESTING NLP SOLVER")
include("NLPSolver/runtests.jl")
println("END TESTING NLP SOLVER")

println("BEGIN TESTING SIP SOLVER")
include("SemiInfinite/runtests.jl")
println("END TESTING SIP SOLVER")
