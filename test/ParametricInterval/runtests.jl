# write your own tests here
println("Testing Checks and Utilities...")
t = @elapsed include("ParamChk_Tests.jl")
println("done (took $t seconds).")

println("Testing Parametric Contractors...")
t = @elapsed include("ParametricContractor_Tests.jl")
println("done (took $t seconds).")
