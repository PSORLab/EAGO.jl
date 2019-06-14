function mc_dense_newton_gs!(z_mc::Vector{MC{N}}, x_mc, YdH_mc::VecOrMat{MC{N}},
                             YH_mc::Vector{MC{N}}, nx::Int) where N


    #println("z_mc: $z_mc")
    #println("x_mc: $x_mc")
    #println("YdH_mc: $YdH_mc")
    #println("YH_mc: $YH_mc")

    S1::MC{N} = zero(x_mc[1])
    x_mc_int::MC{N} = zero(x_mc[1])

    for i=1:nx
      S1 = zero(x_mc[1])
      x_mc_int =  MC{N}(x_mc[i])
      for j=1:nx
        if (i<j)
          S1 = S1 + YdH_mc[i,j]*(x_mc[j] - z_mc[j])
        elseif (j<i)
          S1 = S1 + YdH_mc[i,j]*(x_mc[j] - z_mc[j])
        end
      end
      #println("YH_mc[i] + S1: $(YH_mc[i] + S1)")
      #println("(YH_mc[i] + S1)/YdH_mc[i,i]: $((YH_mc[i] + S1)/YdH_mc[i,i])")
      #println("-(YH_mc[i] + S1)/YdH_mc[i,i]: $(-(YH_mc[i] + S1)/YdH_mc[i,i])")
      eval = z_mc[i] - (YH_mc[i] + S1)/YdH_mc[i,i]
      x_mc[i] = z_mc[i] - (YH_mc[i] + S1)/YdH_mc[i,i]
      #println("eval: $eval")
      x_mc[i] = final_cut(x_mc[i], x_mc_int)
      #eval_cut = final_cut(eval, x_mc_int)
      #println("eval & cut: $eval_cut")
    end

end

function mc_dense_krawczyk_cw!(z_mc::Vector{MC{N}}, x_mc,
                               YdH_mc::VecOrMat{MC{N}}, YH_mc::Vector{MC{N}}, nx::Int) where N

  S1::MC{N} = zero(x_mc[1])
  x_mc_int::MC{N} = zero(x_mc[1])

  for i=1:nx
    S1 = zero(x_mc[1])
    x_mc_int = MC{N}(x_mc[i])
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
