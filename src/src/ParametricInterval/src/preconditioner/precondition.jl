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
        H::Vector{Interval{T}} = h(mid.(X),P)
        J::VecOrMat{Interval{T}} = hj(X,P)
        return H,J
end

"""
    Dense_Precondition!(H,J,Y,opt)
"""
function Dense_Precondition(h::Function,
                            hj::Function,
                            X::Vector{Interval{T}},
                            P::Vector{Interval{T}},
                             opt::PIntvParams) where {T}
    H::Vector{Interval{T}} = h(mid.(X),P)
    J::VecOrMat{Interval{T}} = hj(X,P)
    Y::VecOrMat{T} = mid.(J)
    if (opt.nx == 1)
        YH::Vector{Interval{T}} = H/Y[1,1]
        YJ::VecOrMat{Interval{T}} = J/Y[1,1]
    else
        F = lufact(Y)
        YH = F\H
        YJ = F\J
    end
    return YH,YJ
end

function Dense_Precondition(h::Function,
                            hj::Function,
                            X::Vector{MCInterval{T}},
                            P::Vector{MCInterval{T}},
                            opt::PIntvParams) where {T}
    H::Vector{MCInterval{T}} = h(mid.(X),P)
    J::VecOrMat{MCInterval{T}} = hj(X,P)
    Y::VecOrMat{T} = mid.(J)
    if (opt.nx == 1)
        YH::Vector{MCInterval{T}} = H/Y[1,1]
        YJ::VecOrMat{MCInterval{T}} = J/Y[1,1]
    else
        F = lufact(Y)
        YH = F\H
        YJ = F\J
    end
    return YH,YJ
end
