"""
    param_intv_contractor

Runs an parametic interval contractor method with the options specified in the
opt field. At each iterations, checks based on Miranda are performed for inclusion
and exclusion. Inputs:
    - `h::Function`: where h is the state space function
    - `hj::Function`: where hj is the Jacobian of the state space function w.r.t. x
    - `X::Vector{IntervalType}`: state space variable bounds
    - `P::Vector{IntervalType}`: control space variable bounds
    - `Eflag::Bool`: flag indicating it is known with certainty that a function
                     is exlcuded from the interval box
    - `Iflag::Bool`: flag indicating it is known with certainty that a function
                     is included from the interval box
    - `eDflag::Bool`: flag whether or not extended division occured in the interval box
    - `opt::parametric_interval_params`: option parameters
"""
function param_intv_contractor(h::Function, hj::Function, X, Ntemp::Vector{IntervalType},
                                N::Vector{IntervalType},
                                Xi::Vector{IntervalType}, X1, Xold, t, Y, J, H,
                                P::Vector{IntervalType}, inc::Vector{Bool},
                                incLow::Vector{Bool}, incHigh::Vector{Bool},
                                nx::Int, kmax::Int, etol::Float64, rtol::Float64)

    # Setup storage objects
    fill!(inc, false)
    fill!(incLow, false)
    fill!(incHigh, false)

    k::Int = 1
    eD::Int = 0
    Eflg::Bool = false
    eDflg::Bool = false
    Iflg::Bool = false
    S1 = zero(IntervalType)
    S2 = zero(IntervalType)

    # Runs first iteration
    N[:] = X
    Xi[:] = X
    X1[:] = X
    Eflg, eDflg = parametric_interval_kernel!(h, hj, P, N, Ntemp, X, X1, Xold, t, Y, J, H,
                                              inc, incLow, incHigh, nx, rtol)

    if ~eDflg
        Iflg = inclusion_test(Iflg, inc, nx)
        while ((k < kmax) && (is_equal(X, Xi, etol) == false))
            Xi[:] = X
            Eflg, eDflg = parametric_interval_kernel!(h, hj, P, N, Ntemp, X, X1, Xold,
                                                      t, Y, J, H, inc,
                                                      incLow, incHigh, nx, rtol)
            if (eDflg)
                break
            end
            Iflg = inclusion_test(Iflg, inc, nx)
            k += 1
        end
    end

    if ~eDflg
        X1[:] = X
    end

    return Eflg, Iflg, eDflg
end


"""
    parametric_interval_kernel(X1,X2,S1,S2,H,J,eD,opt)

Contains the kernel of the contractor object. It selects the appropriate one
based on the contractor type selected `CType` and the linear algebra type
selected `LAlg`.
"""
function parametric_interval_kernel!(h::Function, hj::Function,
                                     P::Vector{IntervalType}, N::Vector{IntervalType},
                                     Ntemp::Vector{IntervalType}, X1, X2, Xold,
                                     t::Vector{Float64}, Y, J, H, inc::Vector{Bool},
                                     inclL::Vector{Bool}, inclH::Vector{Bool},
                                     nx::Int, rtol::Float64)
    eDflag::Bool = false
    Eflag::Bool = false
    dense_precondition(h, hj, H, J, X1, Xold, P, t, Y, nx)
    eDflag, Eflag = dense_newton_gs(H, J, N, Ntemp, X1, X2, inc, inclL, inclH, nx, rtol)
    return Eflag, eDflag
end
