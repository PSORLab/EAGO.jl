# write your own tests here
println("Testing Checks and Utilities...")
t = @elapsed include("param_check_tests.jl")
println("done (took $t seconds).")

println("Testing Preconditioners...")
t = @elapsed include("preconditioner_test.jl")
println("done (took $t seconds).")

println("Testing Parametric Contractors...")
t = @elapsed include("parametric_contractor_tests.jl")
println("done (took $t seconds).")
