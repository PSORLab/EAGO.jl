"""
    Precondition!(H,J,Y,opt)

Preconditions the vector H and matrix J by the inverse of the preconditioning
matrix Y.
"""
function Precondition(h::Function,
                      hj::Function,
                      X::Vector{V},
                      P::Vector{V},
                      opt::PIntvParams) where {V<:AbstractInterval}
    if (opt.LAlg == :DenseBand)
        H,J = DenseBand_Precondition!(h,hj,X,P,opt)
    elseif (opt.LAlg == :DenseBlockDiag)
        H,J = DenseBlockDiag_Precondition!(h,hj,X,P,opt)
    elseif (opt.LAlg == :Dense)
        H,J = Dense_Precondition(h,hj,X,P,opt)
    elseif (opt.LAlg == :None)
        H,J = None_Precondition(h,hj,X,P,opt)
    else
        error("Unsupported Linear Algebra Style")
    end
    return H,J
end

function None_Precondition(h::Function,
                           hj::Function,
                           X::Vector{T},
                           P::Vector{T},
                           opt::PIntvParams) where {T<:AbstractInterval}
        H::Vector{T} = h(mid.(X),P)
        J::VecOrMat{T} = hj(X,P)
        return H,J
end

"""
    Dense_Precondition!(H,J,Y,opt)
"""
function Dense_Precondition(h::Function,
                            hj::Function,
                            X::Vector{T},
                            P::Vector{T},
                            opt::PIntvParams{S}) where {T<:AbstractInterval,S<:AbstractFloat}
    println("pre-H")
    H::Vector{T} = h(mid.(X),P)
    println("post H")
    J::VecOrMat{T} = hj(X,P)
    println("post HJ")
    Y::VecOrMat{S} = mid.(J)
    if (opt.nx == 1)
        YH::Vector{T} = H/Y[1,1]
        YJ::VecOrMat{T} = J/Y[1,1]
    else
        F = lufact(Y)
        YH = F\H
        YJ = F\J
    end
    return YH,YJ
end
