# write your own tests here

println("Testing McCormick Multiplication Operator...")
include("Mult_Test.jl")

println("Testing McCormick Operators...")
include("Operators.jl")

println("Testing McCormick Extrema Operator...")
include("Extrema_Test.jl")

println("Testing Utility Functions...")
include("Utilities.jl")

println("Testing McCormick & Constant Operators...")
include("Operators/WithConstants.jl")

println("Implicit Bounding Subroutines...")
include("ImplicitBnd.jl")

#=
println("Implicit Bounding Utilities...")
t = @elapsed include("Imp_Util_Test.jl")
println("done (took $t seconds).")

println("Implicit Bounding Subroutines...")
t = @elapsed include("ImplicitBnd.jl")
println("done (took $t seconds).")
=#

#=
println("Implicit Bounding Utilities...")
t = @elapsed include("D1_Interval_Test.jl")
println("done (took $t seconds).")

println("Implicit Function Routines...")
t = @elapsed include("D1_Interval_Test.jl")
println("done (took $t seconds).")
=#
