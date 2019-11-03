"""
    mc_dense_newton_gs!

Performs a single step of the dense-newton gauss-siedel parametric method after
the inputs have been preconditioned.
"""
function mc_dense_newton_gs!(z_mc::Vector{MC{N,T}}, x_mc, YdH_mc::VecOrMat{MC{N,T}},
                             YH_mc::Vector{MC{N,T}}, nx::Int64) where {N, T<:RelaxTag}

    S1::MC{N,T} = zero(x_mc[1])
    x_mc_int::MC{N,T} = zero(x_mc[1])

    for i=1:nx
      S1 = zero(x_mc[1])
      x_mc_int =  MC{N,T}(x_mc[i])
      for j=1:nx
        if (i<j)
          S1 = S1 + YdH_mc[i,j]*(x_mc[j] - z_mc[j])
        elseif (j<i)
          S1 = S1 + YdH_mc[i,j]*(x_mc[j] - z_mc[j])
        end
      end
      x_mc[i] = z_mc[i] - (YH_mc[i] + S1)/YdH_mc[i,i]
      x_mc[i] = final_cut(x_mc[i], x_mc_int)
    end
end

"""
    mc_dense_krawczyk_cw!

Performs a single step of the dense-krawczyk componentwise parametric method after
the inputs have been preconditioned.
"""
function mc_dense_krawczyk_cw!(z_mc::Vector{MC{N,T}}, x_mc,
                               YdH_mc::VecOrMat{MC{N,T}}, YH_mc::Vector{MC{N,T}},
                               nx::Int64) where {N,T<:RelaxTag}

  S1::MC{N,T} = zero(x_mc[1])
  x_mc_int::MC{N,T} = zero(x_mc[1])

  for i=1:nx
    S1 = zero(x_mc[1])
    x_mc_int = MC{N,T}(x_mc[i])
    for j=1:nx
      if (i<j)
        S1 = S1 - (YdH_mc[i,j])*(x_mc[j] - z_mc[j])
      elseif (j<i)
        S1 = S1 - (YdH_mc[i,j])*(x_mc[j] - z_mc[j])
      else
        S1 = S1 + (one(x_mc[1]) - YdH_mc[i,j])*(x_mc[j] - z_mc[j])
      end
    end
    x_mc[i] =  z_mc[i] - YH_mc[i] + S1
    x_mc[i] = final_cut(x_mc[i],x_mc_int)
  end

end
