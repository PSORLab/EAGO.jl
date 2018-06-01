function SMCg_Dense_Newton_GS!(z_mc::Vector{HybridMC{N,V,T}},x_mc::Vector{HybridMC{N,V,T}},
                        YdH_mc::VecOrMat{HybridMC{N,V,T}},YH_mc::Vector{HybridMC{N,V,T}},
                        mc_opts::mc_opts{T}) where {N,V<:AbstractInterval,T<:AbstractFloat}
    S1::HybridMC{N,V,T} = zero(x_mc[1])
    x_mc_int::Vector{HybridMC{N,V,T}} = copy(x_mc)
    for i=1:mc_opts.nx
      S1 = zero(x_mc[1])
      for j=1:mc_opts.nx
        if (i<j)
          S1 = S1 + YdH_mc[i,j]*(x_mc[j]-z_mc[j])
        elseif (j<i)
          S1 = S1 + YdH_mc[i,j]*(x_mc[j]-z_mc[j])
        end
      end
      x_mc[i] = z_mc[i] - (YH_mc[i]+S1)/YdH_mc[i,i]
      x_mc[i] = Final_Cut(x_mc[i],x_mc_int[i])
    end
end

function SMCg_Dense_Krawczyk_CW!(z_mc::Vector{HybridMC{N,V,T}},x_mc::Vector{HybridMC{N,V,T}},
                                 YdH_mc::VecOrMat{HybridMC{N,V,T}},
                                 YH_mc::Vector{HybridMC{N,V,T}},mc_opts::mc_opts{T}) where{N,V,T<:AbstractFloat}
  S1::HybridMC{N,V,T} = zero(x_mc[1])
  x_mc_int::Vector{HybridMC{N,V,T}} = copy(x_mc)
  for i=1:mc_opts.nx
    S1 = zero(x_mc[1])
    for j=1:mc_opts.nx
      if (i<j)
        S1 = S1 - (YdH_mc[i,j])*(x_mc[j]-z_mc[j])
      elseif (j<i)
        S1 = S1 - (YdH_mc[i,j])*(x_mc[j]-z_mc[j])
      else
        S1 = S1 + (one(T)-YdH_mc[i,j])*(x_mc[j]-z_mc[j])
      end
    end
    x_mc[i] =  z_mc[i] - YH_mc[i] + S1
    x_mc[i] = Final_Cut(x_mc[i],x_mc_int[i])
  end
end
