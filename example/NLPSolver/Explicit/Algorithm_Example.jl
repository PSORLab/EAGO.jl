workspace()

using EAGO
using MathProgBase

# Sets functions defining optimization problem
ffunc(x) = x[1]^3+x[2]+max(x[1],x[2])      # objective function f(x)
ggfunc(x) = [x[1]; -x[2]*x[1]]             # constraint function g(x,p) <= 0 (should input vector and output column vector)
hhfunc(x) = [x[1]-2*x[2]]                  # constraint function h(x,p) == 0 (should input vector and output column vector)

# Sets up solvers which defines multiple problem options
solver_set = EAGO_NLPSolver(LBD_func_relax = "Diff1-MV",      # Use differentiable McCormick relaxations
                            LBDsolvertype = "Ipopt",          # Use Ipopt solver for lower problem
                            UBDsolvertype = "Ipopt",          # Use Ipopt solver for upper problem
                            probe_depth = -1,                 # Disable probing (algorithm under construction)
                            variable_depth = 1000,            # enable DBBT on variable bounds
                            DAG_depth = -1,                   # Disable FW interval constractor (not supported for algorithm interface current, can be used with JuMP front-end)
                            validated = true,                 # Interval calculatons should be correctly rounded (Leave this true for now)
                            atol=1E-4,                        # Sets absolute tolerance
                            rtol=1E-4                         # Sets relative tolerance
                            )

# Sets up and solves constrained explicit problem
# requires bounds on the variables and the objective function, other arguments
# are input as keywords and outputs a EAGO_NLP_Model <: MathProgBase.AbstractNonlinearModel
m = Build_Script(ffunc,                  # Objective function
                 Float64[-10.0,-10.0],   # Lower bounds
                 Float64[10.0,10.0],     # Upper bounds
                 g = ggfunc,             # Inequality constraints
                 h = hhfunc,             # Equality constraints
                 solver = solver_set     # Sets solver
                 )

# If you want to manipulate the internal branch and bound algorithm directly here would be the point
# For instance if you wanted to custom code an function evaluation based upper bound below you could do:
#=
                            """
                                Interval_UBD

                            Solves a upper bounding problem based by overloading functions with correctly
                            rounded natural interval extensions from the ValidatedNumerics.jl packages. Inputs:
                            * `X::Vector{Interval}`: Node over which to solve the upper problem
                            * `k::Int64`: Number of Iterations The B&B Algorithm has taken
                            * `pos::Int64`: Depth of Node in B&B Tree
                            * `opt`: Option type containing problem information
                            * `temp`: The last element of the tuple is currently unused for this option.
                            Returns a tuple `(val,pnt,feas,X,[feas,val])` where
                            * `val::Float64` - Upper bound calculated
                            * `pnt::Vector{Float64}`: An array of length equal to X that gives the optimal
                                                      solution of the upper bound problem.
                            * `feas::Bool`: Returns true if the problem is feasible and false if it is infeasible
                            """
                            function Interval_UBD(X::Vector{V},k::Int,pos::Int,opt,temp) where {V<:AbstractInterval}

                                  # solve optimization problem via interval extension
                                  feas = true
                                  pnt::Vector{Float64} = mid.(X)
                                  FInt::V = opt[1].f(V.(pnt))
                                  val::Float64 = FInt.hi
                                  (opt[1].numConstr > 0) && error("Problem must be unconstrained for interval-midpoint upper bound.")

                                  return val, pnt, feas, Any[feas,val]
                            end

                            m.Opts.solver.Upper_Prob = Interval_UBD
=#


MathProgBase.optimize!(m)  # Optimizes the problem

# extracts objective value, solution, feasibility
stat = MathProgBase.status(m)
soln = MathProgBase.getsolution(m)
objv = MathProgBase.getobjval(m)
