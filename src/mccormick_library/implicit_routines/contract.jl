"""
$(FUNCTIONNAME)

Performs a single step of the parametric method associated with `t` the inputs have been preconditioned.
"""
function contract! end
function contract!(t::NewtonGS, d::MCCallback{FH,FJ,C,PRE,N,T}) where {FH <: Function,
                                                                       FJ <: Function,
                                                                       C, PRE, N,
                                                                       T<:RelaxTag}
    S = zero(MC{N,T})
    @. d.x0_mc = d.x_mc
    for i = 1:d.nx
        S = zero(MC{N,T})
        for j = 1:d.nx
            if (i !== j)
                @inbounds S += d.J[i,j]*(d.x_mc[j] - d.z_mc[j])
            end
        end
        @inbounds d.x_mc[i] = d.z_mc[i] - (d.H[i] + S)*McCormick.inv1(d.J[i,i], 1.0/d.J[i,i].Intv)
        @inbounds d.x_mc[i] = final_cut(d.x_mc[i], d.x0_mc[i])
    end
    return
end
# contract is optimized!
function contract!(t::KrawczykCW, d::MCCallback{FH,FJ,C,PRE,N,T}) where {FH <: Function,
                                                                         FJ <: Function,
                                                                         C, PRE, N,
                                                                         T<:RelaxTag}
    S1::MC{N,T} = zero(MC{N,T})
    x_mc_int::MC{N,T} = zero(MC{N,T})
    for i=1:nx
        S1 = zero(x_mc[1])
        x_mc_int = MC{N,T}(x_mc[i])
        @inbounds @simd for j=1:nx
            if (i<j)
                S1 -= (J[i,j])*(x_mc[j] - z_mc[j])
            elseif (j<i)
                S1 -= (J[i,j])*(x_mc[j] - z_mc[j])
            else
                S1 += (one(x_mc[1]) - J[i,j])*(x_mc[j] - z_mc[j])
            end
        end
        x_mc[i] =  z_mc[i] - H[i] + S1
        x_mc[i] = final_cut(x_mc[i],x_mc_int)
    end
    return
end
# contract is optimized!
