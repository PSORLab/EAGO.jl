#=
"""
    Precondition!(H,J,Y,opt)

Preconditions the vector H and matrix J by the inverse of the preconditioning
matrix Y.
"""
function precondition(h::Function, hj::Function, X::Vector{IntervalType},
                      P::Vector{IntervalType}, opt::parametric_interval_params)
        H,J = dense_precondition(h,hj,X,P,opt)
    elseif (opt.LAlg == :None)
        H,J = none_precondition(h,hj,X,P,opt)
    else
        error("Unsupported Linear Algebra Style")
    end
    return H,J
end
=#

function none_precondition(h!::Function, hj!::Function, X::Vector{IntervalType},
                           P::Vector{IntervalType}, opt::parametric_interval_params)
        H = zeros(IntervalType, opt.nx)
        J = zeros(IntervalType, opt.nx, opt.nx)
        h!(H, mid.(X), P)
        hj!(J, X, P)
        return H,J
end

"""
    dense_precondition!(H,J,Y,opt)
"""
function dense_precondition(h!::Function, hj!::Function, H, J, X, Xold, P::Vector{IntervalType},
                            t::Vector{Float64}, Y::Array{Float64,2}, nx::Int)

    fill!(H, zero(IntervalType))
    fill!(J, zero(IntervalType))

    xmid = mid.(X)

    h!(H, xmid, Xold, P, t)
    hj!(J, X, Xold, P, t)


    Y .= mid.(J)
    if (nx == 1)
        H[1] = H[1]/Y[1,1]
        J[1,1] = J[1,1]/Y[1,1]
    else
        #Yinv = inv(Y)
        #H[:] = Yinv*H
        #J[:,:] = Yinv*J
        Y[:,:] = inv(Y)
        H[:] = Y*H
        J[:,:] = Y*J
    end
end
