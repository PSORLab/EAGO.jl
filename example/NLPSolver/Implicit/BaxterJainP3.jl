
function Pressure_Matrix(Lpt::T,Kt::T,N::Int,R::Float64,Svt::Float64,Pvv::Float64)

    # Pressure Profile
    att = R*sqrt(Lpt*Svt/Kt)

    # N -- number of spatial grid points
    r = linspace(0,R,N)
    r = r/R
    dr = 1/(N-1)
    ivdr2 = 1./dr^2
    att2 = att^2

    A = zeros(N,N)

    for i in 2:N-1
        A[i,i-1] = -1/r[i]/dr + ivdr2
        #A[i,i] = -2./dr^2 - (at[i]^2)
        A[i,i] = -2.*ivdr2 - (att2)
        A[i,i+1] = 1./r[i]/dr + ivdr2
        #F[i] = -(at[i]^2)*Pvv
        F[i] = -(att2)*Pvv
    end

    # Boundary condtions for isolated tumor model
    A[1,1] = -2/dr^2 - (att^2)
    A[1,2] = 2/dr^2
    A[N,N] = 1.0
end

function Pressure_RHS(Lpt::T,Kt::T,N::Int,R::Float64,Svt::Float64,Pvv::Float64)

    # Pressure Profile
    att = R*sqrt(Lpt*Svt/Kt)

    # N -- number of spatial grid points
    r = linspace(0,R,N)
    r = r/R
    dr = 1/(N-1)
    ivdr2 = 1./dr^2
    att2 = att^2

    F = zeros(N,1)

    for i in 2:N-1
        F[i] = -(att2)*Pvv
    end

    F[1] = -att2*Pvv
end

function Concentration_Matrix(P,V,)

end

#
function Conc_RHS(i::Int,P::,h)

    co::Float64 = 1.0                                 # dimensionless drug concentration
    tspan::Float64 = 1.0*3600.0                       # length of simulation type
    cv::Float64 = co*exp(-t/kd/3600.0)                # vascular concentration of the drug following exponential decay
    coeff2::T = h*Lpt*(Svt*Pv*cv*(1.0-sigma))
    coeff1::Float64 = h*Peff*Svt*cv+coeff2*Pvv

    f::Vector{T} = coeff1-coeff2*P
    f[end] = zero(T)

    return f
end

function LowerBound(X::Vector{Interval{Float64}},k::Int64,pos::Int64,opt::Any,UBD::Float64)

    opt[1].solver.SubGradRefine && set_hybrid_box!(X,x0,true)
    
    N = 100
    R = 1.0
    Svt = 1 # fix me
    Pvv = 1 # fix me
    opt = mc_opt(Float64)

    p = 1
    pref = 1

    # Get bounds for Kt and Lpt

    # Bound pressure via implicit form
    Pressure_M = Pressure_Matrix(Lpt,Kt,N,R,Svt,Pvv)
    Pressure_RHS = Pressure_RHS(Lpt,Kt,N,R,Svt,Pvv)
    P = BandedImplicitGS(p,pref,X,Pressure_M,Pressure_RHS,opt)

    # Calculate V bounds
    V = vcat(zero(T),diff(P))

    # Bound concentration via implicit form
    Conc_M = Concentration_Matrix(P,V,)
    Conc_M = inv(Conc_M)*Conc_M
    rhs = Conc_RHS(i+1)
    for i=1:99
        Conc_RHS = rhs + X[(N*(i-1)+1):(N*(i-1)+N)]
        X[(N*i+1:(N*i+N)] = BandedImplicitGS(p,pref,X[(N*i+1:(N*i+N)],Conc_M,Conc_RHS,opt)
    end
end
