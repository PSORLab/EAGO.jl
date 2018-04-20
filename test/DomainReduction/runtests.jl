println("Testing Duality-Based Bound Tightening...")
t = @elapsed include("TestDual.jl")
println("done (took $t seconds).")

println("Testing Range Reduction...")
t = @elapsed include("TestRR.jl")
println("done (took $t seconds).")

println("Testing DAG Contractor ...")
t = @elapsed include("TestDAGcntr.jl")
println("done (took $t seconds).")
