"""
"""
function Dense_Newton_GS(H::Vector{Th},J::VecOrMat{Tj},N::Vector{V},
                         S1::V,S2::V,X1::Vector{V},X2::Vector{V},
                         incl::Vector{Bool},inclL::Vector{Bool},
                         inclH::Vector{Bool},opt::PIntvParams{T}) where {Th,Tj,V,T}
    for i=1:opt.nx
        S1 = zero(V)
        S2 = zero(V)
        for j=1:opt.nx
            if (j<i)
                S1 += J[i,j]*(X1[j]-mid.(X1[j]))
            elseif (j>i)
                S2 += J[i,j]*(X1[j]-mid.(X1[j]))
            end
        end
        #println("i: $i")
        #println("check val: $(J[i,i].lo*J[i,i].hi)")
        if J[i,i].lo*J[i,i].hi > 0.0
          N[i] = mid(X1[i]) - (H[i]+S1+S2)/J[i,i]
          #println("X1[i]: $(X1[i])")
          #println("N[i]: $(N[i])")
        else
          Ntemp = copy(N)
          eD,N[i],Ntemp[i] = extProcess(N[i],X1[i],J[i,i],S1,S2,H[i],opt.rtol)
         # println("eD: $eD")
         # println("N[i]: $(N[i])")
         # println("Ntemp[i]: $(Ntemp[i])")
          if eD == 1
            eDflg = true
            X2 = copy(X1)
            X2[i] = Ntemp[i] ∩ X1[i]
            X1[i] = N[i] ∩ X1[i]
            return true,false
          end
        end
        Strictly_InRoutine!(i,N,X1,incl,inclL,inclH)
        if ~isdisjoint(N[i],X1[i])
            X1[i] = N[i] ∩ X1[i]
        else
            return false,true
        end
    end
    return false,false
end

#=
"""
"""
function DenseBlockDiag_Newton_GS(H::Vector{S},J::Array{S,2},N::V,
                               S::V,X1::Vector{V},X2::Vector{V},
                               incl::Vector{Bool},inclL::Vector{Bool},
                               inclH::Vector{Bool},opt::PIntvParams{T})
    blk_cnt::Integer = 0
    intr_cnt::Integer = 1
    for i=1:opt.nx
        S1 = zero(V)
        S2 = zero(V)
        for j=(blk_cnt*opt.tblock+1):((blk_cnt+1)*opt.tblock)
            if (j<i)
                S1 += J.data[intr_cnt,j]*(X1[j]-mid.(X1[j]))
            elseif (j>i)
                S2 += (J.data[intr_cnt,j]*(X1[j]-mid.(X1[j]))
            end
        end
        midid = 1+(opt.tblock+1)*(i-1)
        if J.data[midid].lo*J.data[midid].hi > 0.0
          N[i] = mid(X1[i]) - (H[i]+S1+S2)/J.data[midid]
        else
          Ntemp = copy(N)
          eD,N[i],Ntemp[i] = extProcess(N[i],X1[i],J.data[midid],S1,S2,H[i],opt.rtol)
          if eD == 1
            eDflg = true
            X2 = copy(X1)
            X2[i] = Ntemp[i] ∩ X1[i]
            X1[i] = N[i] ∩ X1[i]
            return true,false
          end
        end
        Strictly_InRoutine!(i,N,X1,incl,inclL,inclH)
        if ~isdisjoint(N[i],X1[i])
            X1[i] = N[i] ∩ X1[i]
        else
            return false,true
        end
        intr_cnt += 1
        if (intr_cnt > opt.tblock)
            intr_cnt = 1
            blk_cnt += 1
        end
    end
    return false,false
end

"""
"""
function DenseBand_Newton_GS(H::Vector{S},J::Array{S,2},N::V,
                               S::V,X1::Vector{V},X2::Vector{V},
                               incl::Vector{Bool},inclL::Vector{Bool},
                               inclH::Vector{Bool},opt::PIntvParams{T})

    blk_cnt::Integer = 0
    intr_cnt::Integer = 1
    for i=1:opt.nx
        S1 = zero(V)
        S2 = zero(V)
        for j=(blk_cnt*opt.tblock+1):((blk_cnt+1)*opt.tblock)
            if (j<i)
                S1 += J.data[intr_cnt,j]*(X1[j]-mid.(X1[j]))
            elseif (j>i)
                S2 += (J.data[intr_cnt,j]*(X1[j]-mid.(X1[j]))
            end
        end
        midid = 1+(opt.tblock+1)*(i-1)
        if J.data[midid].lo*J.data[midid].hi > 0.0
            N[i] = mid(X1[i]) - (H[i]+S1+S2)/J.data[midid]
        else
            Ntemp = copy(N)
            eD,N[i],Ntemp[i] = extProcess(N[i],X1[i],J.data[midid],S1,S2,H[i],opt.rtol)
            if eD == 1
                eDflg = true
                X2 = copy(X1)
                X2[i] = Ntemp[i] ∩ X1[i]
                X1[i] = N[i] ∩ X1[i]
                return true,false
             end
         end
         Strictly_InRoutine!(i,N,X1,incl,inclL,inclH)
         if ~isdisjoint(N[i],X1[i])
             X1[i] = N[i] ∩ X1[i]
         else
             return false,true
         end
         intr_cnt += 1
         if (intr_cnt > opt.tblock)
             intr_cnt = 1
             blk_cnt += 1
         end
    end
    return false,false
end
=#
