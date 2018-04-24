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
