# Test basic interval arithmetic functions
println("Testing Arithmetic...")
t = @elapsed include("arithmetic_test.jl")
println("done (took $t seconds).")
