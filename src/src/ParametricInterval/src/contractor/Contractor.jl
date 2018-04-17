"""
    Param_Intv_Contractor

Runs an parametic interval contractor method with the options specified in the
opt field. At each iterations, checks based on Miranda are performed for inclusion
and exclusion.

Returns
    - `h::Function`: where h is the state space function
    - `hj::Function`: where hj is the Jacobian of the state space function
                      w.r.t. x
    - `X::Vector{S}`: where S is either Interval or MCInterval
    - `P::Vector{S}`: where S is either Interval or MCInterval
    - `Eflag::Bool`: flag indicating it is known with certainty that a function
                     is exlcuded from the interval box
    - `Iflag::Bool`: flag indicating it is known with certainty that a function
                     is included from the interval box
    - `eDflag::Bool`: flag whether or not extended division occured in the
                      interval box
    - `opt::PIntvParams{T}`: option parameters
"""
function Param_Intv_Contractor(h::Function,hj::Function,
                               X::Vector{S},P::Vector{S},
                               Eflag::Bool,Iflag::Bool,
                               eDflag::Bool,opt::PIntvParams{T}) where {S,T}

    # Setup storage objects
    inc::Vector{Bool} = Bool[false for i=1:opt.nx]       # Total inclusion test storage
    incLow::Vector{Bool} = Bool[false for i=1:opt.nx]    # Lower inclusion test storage
    incHigh::Vector{Bool} = Bool[false for i=1:opt.nx]   # Upper inclusion test storage
    k::Int64 = 1                                         # Iteration number
    eD::Int64 = 0                                        # Extended division flag
    Eflg::Bool = false                                   # The exclusion of the function is not known with certainty
    eDflg::Bool = false                                  # The inclusion of the function is not known with certainty
    Iflg::Bool = false                                   # Division by zero occured and extended division was necessary
    S1 = zero(S)
    S2 = zero(S)

    # Runs first iteration
    N::Vector{S} = copy(X)
    Xi::Vector{S} = copy(X)
    X1::Vector{S} = copy(X)
    #println("X start: $X")
    X,X1,Eflg,eDflg,inc = PIntv_Kernel(h,hj,P,N,X,X1,S1,S2,eD,inc,incLow,incHigh,opt)  # performs the contractor
    #println("X end: $X")
    #println("inc at k = $k: $inc")
    (eDflag) && (return X,X1,Eflg,Iflg,eDflg,incLow,incHigh)                      # returns function if extended division occurs
    Iflg = InclTest(Iflg,inc,opt.nx)
    #println("Iflg: $Iflg")                                              # checks for inclusion

    while ((k<opt.kmax) && (isEqual(X,Xi,opt.etol) == false))
        # Precondition H and J using interval midpoint inverse Jacobian
        Xi = copy(X)
        #println("X start: $X")
        X,X1,Eflg,eDflg,inc = PIntv_Kernel(h,hj,P,N,X,X1,S1,S2,eD,inc,incLow,incHigh,opt)
        #println("X end: $X")
        #println("inc at k = $k: $inc")
        (eDflag) && (return X,X1,Eflg,Iflg,eDflg,incLow,incHigh)     # returns function if extended division occurs
        Iflg = InclTest(Iflg,inc,opt.nx)                             # checks for inclusion
        #println("Iflg: $Iflg")
        k+=1
    end

    X1 = copy(X)    # make a copy of X since no extended division occured
    return X,X1,Eflg,Iflg,eDflg,incLow,incHigh
end


"""
    PIntv_Kernel(X1,X2,S1,S2,H,J,eD,opt)

Contains the kernel of the contractor object. It selects the appropriate one
based on the contractor type selected `CType` and the linear algebra type
selected `LAlg`.
"""
function PIntv_Kernel(h::Function,hj::Function,
                      P::Vector{S},N::Vector{S},
                      X1::Vector{S},X2::Vector{S},
                      S1::S,S2::S,eD::Integer,
                      incl::Vector{Bool},inclL::Vector{Bool},inclH::Vector{Bool},
                      opt::PIntvParams) where {S}
    eDflag::Bool = false
    Eflag::Bool = false
    H,J = Precondition(h,hj,X1,P,opt)
    if (opt.CTyp == :Newton)
        if (opt.LAlg == :Dense)
            eDflag,Eflag = Dense_Newton_GS(H,J,N,S1,S2,X1,X2,incl,inclL,inclH,opt)
        elseif (opt.LAlg == :DenseBand)
            eDflag,Eflag = DenseBand_Newton_GS(H,J,N,S1,S2,X1,X2,incl,inclL,inclH,opt)
        else
            error("The linear algebra type $(LAlg) is not currently supported. The
                   linear algebra styles currently supported are :Dense and
                   :DenseBanded.")
        end
    elseif (opt.CTyp == :Krawczyk)
        if (opt.LAlg == :Dense)
            eDflag,Eflag = Dense_Krawczyk_CW(H,J,N,S1,X1,X2,incl,inclL,inclH,opt)
        elseif (opt.LAlg == :DenseBand)
            eDflag,Eflag = DenseBand_Krawczyk_CW(H,J,N,S1,S2,X1,X2,incl,inclL,inclH,opt)
        else
            error("The linear algebra type $(LAlg) is not currently supported. The
                   linear algebra styles currently supported are :Dense and
                   :DenseBanded.")
        end
    else
        error("The contractor type $(CType) is not currently supported. The
               contractors :Newton and :Krawczyk are currently supported.")
    end

    return X1,X2,Eflag,eDflag,incl
end
