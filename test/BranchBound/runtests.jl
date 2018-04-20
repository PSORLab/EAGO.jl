
# Runs a simple optimization problem using Interval Bounds
println("Testing 1D Interval Optimization...")
t = @elapsed include("D1_Interval_Test.jl")
println("done (took $t seconds).")

# Sets options
println("Testing Option Setting...")
t = @elapsed include("Options_Tests.jl")
println("done (took $t seconds).")

# Test node selection routines
println("Testing Nodes...")
t = @elapsed include("Node_Tests.jl")
println("done (took $t seconds).")

# Test bisection routines
println("Testing Bisection...")
t = @elapsed include("Bisect_Tests.jl")
println("done (took $t seconds).")

# Test branching routine
println("Testing Branch...")
t = @elapsed include("Branch_Tests.jl")
println("done (took $t seconds).")

# Tests access routines
println("Testing Access...")
t = @elapsed include("Access_Tests.jl")
println("done (took $t seconds).")

# Tests various checks
println("Testing Checks...")
t = @elapsed include("Checks_Tests.jl")
println("done (took $t seconds).")

# Tests displaying information
println("Testing Display...")
t = @elapsed include("Display_Tests.jl")
println("done (took $t seconds).")
