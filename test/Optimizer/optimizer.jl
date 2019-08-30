println("BEGIN TESTING MOI OPTIMIZER INTERFACE...")
#include("moi_interface.jl")
println("TESTING  MOI OPTIMIZER INTERFACE COMPLETE.")

println("BEGIN TESTING QUADRATIC RELAXATIONS...")
#include("quadratic_relaxation.jl")
println("TESTING QUADRATIC RELAXATIONS COMPLETE.")

println("BEGIN TESTING STANDARD EVALUATOR...")
#include("standard_evaluator.jl")
println("TESTING STANDARD EVALUATOR COMPLETE.")

println("BEGIN TESTING IMPLICIT EVALUATOR...")
#include("implicit_optimizer.jl")
println("TESTING IMPLICIT EVALUATOR COMPLETE.")

println("BEGIN LP TEST PROBLEMS...")
for i in 1:4
   include("TestProblems/LP/Prob$i.jl")
end
println("LP TEST COMPLETE.")

println("BEGIN QP TEST PROBLEMS...")
for i in 1:1
#        include("TestProblems/QP/Prob$i.jl")  # ADD 2 and 3
end
println("QP TEST COMPLETE.")


println("BEGIN NLP TEST PROBLEMS...")
for i in 1:5
   include("TestProblems/NLP/Prob$i.jl")
end
println("NLP TEST COMPLETE.")
