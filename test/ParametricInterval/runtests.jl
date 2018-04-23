# write your own tests here
println("Testing Checks and Utilities...")
t = @elapsed include("ParamChk_Tests.jl")
println("done (took $t seconds).")

println("Testing Preconditioners...")
t = @elapsed include("Preconditioner_Test.jl")
println("done (took $t seconds).")

println("Testing Parametric Contractors...")
t = @elapsed include("ParametricContractor_Tests.jl")
println("done (took $t seconds).")
