# C:\Users\wilhe\AppData\Local\Programs\Julia-1.6.0\bin\julia.exe --track-allocation=user 
#

using JuMP, EAGO

m = Model(EAGO.Optimizer)
set_optimizer_attribute(m, "output_iterations", 1)
set_optimizer_attribute(m, "iteration_limit", 5)
set_optimizer_attribute(m, "verbosity", 0)

# Define bounded variables
xL = [10.0; 0.0; 0.0; 0.0; 0.0; 85.0; 90.0; 3.0; 1.2; 145.0]
xU = [2000.0; 16000.0; 120.0; 5000.0; 2000.0; 93.0; 95.0; 12.0; 4.0; 162.0]
@variable(m, xL[1] <= x[1] <= xU[1])
@NLobjective(m, Max, -2.0 * x[1])


JuMP.optimize!(m)