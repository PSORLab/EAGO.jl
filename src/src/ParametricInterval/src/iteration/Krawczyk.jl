"""
    Dense_Krawczyk_CW

Applies a Krawczyk iteration componentwise using a dense data format for J.
"""
function Dense_Krawczyk_CW(H::Vector{Th},J::Array{Tj,2},N::Vector{V},
                           S::V,X1::Vector{V},X2::Vector{V},
                           incl::Vector{Bool},inclL::Vector{Bool},
                           inclH::Vector{Bool},opt::PIntvParams{T}) where {Tj,Th,V,T}
    for i=1:opt.nx
        S = zero(V)
        for j=1:opt.nx
            if (j != i)
                S += -J[i,j]*(X1[j]-mid.(X1[j]))
            elseif (j == i)
                S += (one(T)-J[i,j])*(X1[j]-mid.(X1[j]))
            end
        end
        N[i] = mid(X1[i]) - H[i] + S
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
    DenseBand_Krawczyk_CW

Mutates inclL, inclH, X1, X2 in place.
Returns an extended division (always false) and exclusion flag.
"""
function DenseBand_Krawczyk_CW(H::Vector{S},J::Array{S,2},N::V,
                               S::V,X1::Vector{V},X2::Vector{V},
                               incl::Vector{Bool},inclL::Vector{Bool},
                               inclH::Vector{Bool},opt::PIntvParams{T}) where {S,V,T}
    blk_cnt::Integer = 0
    intr_cnt::Integer = 1
    for i=1:opt.nx
        S = zero(V)
        for j=(blk_cnt*opt.tblock+1):((blk_cnt+1)*opt.tblock)
            if (j<i)
                S -= J.data[intr_cnt,j]*(X1[j]-mid.(X1[j]))
            elseif (j==i)
                S += (one(T)-J.data[intr_cnt,j])*(X1[j]-mid.(X1[j]))
            end
        end
        N[i] = mid(X1[i]) - H[i] + S
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

#=
"""
    DenseBlockDiag_Krawczyk_CW

Mutates inclL, inclH, X1, X2 in place.
Returns an extended division (always false) and exclusion flag.
"""
function DenseBlockDiag_Krawczyk_CW(H::Vector{S},J::Array{S,2},N::V,
                               S::V,X1::Vector{V},X2::Vector{V},
                               incl::Vector{Bool},inclL::Vector{Bool},
                               inclH::Vector{Bool},opt::PIntvParams{T}) where {S,V,T}
    blk_cnt::Integer = 0
    intr_cnt::Integer = 1
    for i=1:opt.nx
        S = zero(V)
        for j=(blk_cnt*opt.tblock+1):((blk_cnt+1)*opt.tblock)
            if (j<i)
                S -= J.data[intr_cnt,j]*(X1[j]-mid.(X1[j]))
            elseif (j==i)
                S += (one(T)-J.data[intr_cnt,j])*(X1[j]-mid.(X1[j]))
            end
        end
        N[i] = mid(X1[i]) - H[i] + S
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
