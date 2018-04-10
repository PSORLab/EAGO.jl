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
function Interval_UBD(X::Vector{Interval{Float64}},k::Int64,pos::Int64,opt,temp)

      # solve optimization problem via interval extension
      feas = true
      FInt::Interval = opt[1].f(X)
      if (opt[1].numConstr < 1)
      else
          GInt::Vector{Interval} = opt[1].g(X)
          cInt::Vector{Interval{Float64}} = vcat(GInt[opt[1].gU_loc]-opt[1].gU[opt[1].gU_loc],-GInt[opt[1].gL_loc]+opt[1].gL[opt[1].gL_loc])
          if (opt[1].numConstr < 1)
          else
            for i=1:opt[1].numConstr
                if (cInt[i].lo>0.0)
                    feas = false
                    break
                end
            end
          end
      end
      val::Float64 = FInt.hi
      pnt::Vector{Float64} = mid.(X)

       # output
      return val, pnt, feas, Any[feas,val]
end
