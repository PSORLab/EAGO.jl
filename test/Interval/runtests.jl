# Test basic interval arithmetic functions
println("Testing Arithmetic...")
t = @elapsed include("arithmetic_test.jl")
println("done (took $t seconds).")

println("Testing Trig Operations...")
t = @elapsed include("trignometric_test.jl")
println("done (took $t seconds).")

println("Testing Hyperbolic Operations...")
t = @elapsed include("hyperbolic_test.jl")
println("done (took $t seconds).")

println("Testing Set Operations...")
t = @elapsed include("setops.jl")
println("done (took $t seconds).")
