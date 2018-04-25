"""
    Interval_LBD

Solves a lower bounding problem based by overloading functions with correctly
rounded natural interval extensions from the ValidatedNumerics.jl packages.
Inputs to the function are below:
* `X::Vector{Interval{Float64}}`: Node over which to solve the lower problem
* `k::Int64`: Number of Iterations The B&B Algorithm has taken
* `pos::Int64`: Depth of Node in B&B Tree
* `opt`: Option type containing problem information
* `UBD::Float64`: Global upper bound for B&B algorithm
Returns a tuple (val,pnt,feas,X,[]) where
* `val::Float64`: Lower bound calculated
* `pnt::Array{Float64,1}`: An array of length equal to X that gives the
                           optimal solution of the lower bound problem.
* `feas::Bool`: Returns true if the problem is feasible and false if it is
                infeasible
* `[]`:        The last element of the tuple is currently unused for this option.
"""
function Interval_LBD(X::Vector{Q},k::Int64,pos::Int64,opt,UBD) where {Q<:AbstractInterval}

      # solve optimization problem via interval extension
      FInt::Q = opt[1].f(X)
      feas::Bool = true
      if (opt[1].numConstr < 1)
      else
        GInt::Vector{Q} = opt[1].g(X)
        cInt::Vector{Q} = vcat(GInt[opt[1].gU_loc]-opt[1].gU[opt[1].gU_loc],
                                              -GInt[opt[1].gL_loc]+opt[1].gL[opt[1].gL_loc])
        for i=1:opt[1].gexp
          if (cInt[i].lo>0.0)
            feas = false
            break
          end
        end
      end
      val::Float64 = FInt.lo
      pnt::Vector{Float64} = mid.(X)

  # output
  return val, pnt, feas, []
end
