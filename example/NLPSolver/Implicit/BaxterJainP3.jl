
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

function Concentration_Matrix(P,V,Peff,sigma,Svt,h,kd,D,N)

    dr = 1/(N-1)
    Pv::Float64 = 25.0
    co::Float64 = 1.0 # dimensionless drug concentration
    coeff1::T = Peff*Svt
    cv::Float64 = co*exp(-t/kd/3600.0)
    coeff2::T = Lpt*(Svt*Pv*cv*(1.0-sigma))
    coeff3 = 2.0*D/dr
    coeff4 = D/dr^2
    M = zeros(T,N*N)
    M[1,1] = 1.0 + h*(2.0*D/dr^2 + Peff*Svt)
    M[1,2] = -h*2.0*D/dr^2
    for i=2:N
        M[i,i-1] = 1 - h*(coeff4*(c[i-1]))
        M[i,i] = 1 - h*(-(V[i]/dr)-coeff1-2.0*coeff4-(coeff3/r[i]))
        M[i,i+1] = 1 - h*(coeff3/r[i]+coeff4+(V[i]/dr))
    end
    return M
end

#
function Conc_RHS(i::Int,Peff,sigma,Svt,h,kd,D,N)

    Pvv::Float64 = 1.0
    Pv::Float64 = 25.0
    co::Float64 = 1.0                                 # dimensionless drug concentration
    t = i*h
    cv::Float64 = co*exp(-t/kd/3600.0)                # vascular concentration of the drug following exponential decay
    coeff2::T = h*Lpt*(Svt*Pv*cv*(1.0-sigma))
    coeff1::Float64 = h*Peff*Svt*cv+coeff2*Pvv

    f::Vector{T} = coeff1-coeff2*P
    f[end] = zero(T)

    return f
end

function solutePerm(Lpt,rs)

    # calculate diffusion coefficient from Stoke's Einstein
    kB = 1.380648*10.0^(-23)               # Boltzmann Constant (J/K)
    T = 310.15                           # Temperature K
    eta = 3*10.0^(-5)                      # viscosity of blood (mmHg-sec)
    conv_mu  = 133.322365                # (Pascal/mmHg)
    etac = eta*conv_mu                   # Pascal-sec
    pore_conv = 10.0^(-9)                  # (m/nm)
    r_partc = rs*pore_conv               # radius (m)
    D0 = kB*T/(6*pi*etac*r_partc)*1.0e4   # Diffusivity (cm^2/s)

    # Bungay and Brenner
    a = [-73/60,77293/50400,-22.5083,-5.6117,-0.3363,-1.216,1.647]
    b = [7/60;-2227/50400;4.0180;-3.9788;-1.9215;4.392;5.006]

    # Calculate the pore size
    gamma = 1.0e-3
    eta = 3.0e-5                              # Blood viscosity (mmHg/s)
    L = 5.0e-4                                # Vessel wall thickness (cm)
    r_pore = sqrt(8*eta*L*Lpt/gamma)*1.0e7    # nm
    # rs = 30;                              # solute radius (nm) 30nm for FITC
    lambda = rs/r_pore
    t1 = zero(rs)
    t2 = zero(rs)
    p1 = zero(rs)
    p2 = zero(rs)
    for i = 1:7
        if i<3
            t1 = t1 + a[i]*(1-lambda)^i
            p1 = p1 + b[i]*(1-lambda)^i
        else
            t2 = t2 + a[i]*lambda^(i-3)
            p2 = p2 + b[i]*lambda^(i-3)
        end
    end
    Kt = t2 + 9.0/4.0*pi^2*sqrt(2.0)*(1.0+t1)*sqrt(1.0-lambda)^(-5)
    Ks = p2 + 9.0/4.0*pi^2*sqrt(2.0)*(1.0+p1)*sqrt(1.0-lambda)^(-5)
    Phi = (1.0-lambda)^2
    H = 6.0*pi*Phi/Kt
    W = Phi*(2.0-Phi)*Ks/(2.0*Kt)
    Perm = gamma*H*D0/L
    sigma = 1.0 - W
    return Perm,sigma
end

function SMCg_BlockBand_PI_Newton_GS!(z_mc::Vector{HybridMC{N,V,T}},x_mc::Vector{HybridMC{N,V,T}},
                                   YdH_mc::Vector{VecOrMat{HybridMC{N,V,T}}},YH_mc::Vector{Vector{HybridMC{N,V,T}}},
                                   mc_opts::mc_opts{T}) where {N,V<:AbstractInterval,T<:AbstractFloat}
    S1::HybridMC{N,V,T} = zero(x_mc[1])
    row::Int = 0
    x_mc_int::Vector{HybridMC{N,V,T}} = copy(x_mc)
    block_size::Int = size(YdH_mc[1],2)
    block_num::Int = (mc_opts.nx/block_num)
    for k=1:block_num
        for q=1:block_size
            row = (k-1)*block_size + q
            S1 = zero(x_mc[1])
            for j=1:block_size
              if (q<j)
                S1 = S1 + YdH_mc[k][q,j]*(x_mc[row]-z_mc[row])
            elseif (j<q)
                S1 = S1 + YdH_mc[k][q,j]*(x_mc[row]-z_mc[row])
              end
            end
            x_mc[row] = z_mc[row] - (YH_mc[k][row]+S1)/YdH_mc[k][row,row]
            x_mc[row] = Final_Cut(x_mc[row],x_mc_int[row])
        end
    end
end

#=
function SMCg_DenseBlockDiag_Precondition!(h,hj,z_mc,aff_mc,p_mc,opt)
    #Y =
    #hind =
    H = h(z_mc,p_mc)
    J = hj(aff_mc,p_mc)
    block_size::Int = size(J[1],2)
    block_num::Int = (mc_opts.nx/block_num)
    for k=1:block_num
        Y = inv([mid(J[k][i,j].SMC.Intv) for i=1:opt.size, j=1:opt.size])
        J[k] = invY*J[k]
        hind = ((k-1)*block_size+1):(block_size*k)
        H[hind] = invY*H[hind]
    end
    return H,J
end
=#

function Generate_Affine_BlockBand_PI(x_mc::Vector{HybridMC{N,V,T}},
                                      a_mc::VecOrMat{HybridMC{N,V,T}},
                                      b_mc::Vector{HybridMC{N,V,T}},
                                      pmid::Vector{T},mc_opts::mc_opts{T}) where {N,T<:AbstractFloat,V<:AbstractInterval}

  nxi::Int64 = length(X)
  np::Int64 = length(P)

  x_mc::Vector{HybridMC{np,V,T}} = HybridMC{np,V,T}[HybridMC{np,V,T}(SMCg{np,V,T}(X[i].hi,X[i].lo,szero,szero,V(X[i].lo,X[i].hi),false)) for i=1:nxi]
  xa_mc::Vector{HybridMC{np,V,T}} = HybridMC{np,V,T}[HybridMC{np,V,T}(SMCg{np,V,T}(X[i].lo,X[i].lo,szero,szero,V(X[i].lo,X[i].lo),false)) for i=1:nxi]
  xA_mc::Vector{HybridMC{np,V,T}} = HybridMC{np,V,T}[HybridMC{np,V,T}(SMCg{np,V,T}(X[i].hi,X[i].hi,szero,szero,V(X[i].hi,X[i].hi),false)) for i=1:nxi]
  z_mct::Vector{HybridMC{np,V,T}} = mc_opts.lambda*xa_mc+(one(T)-mc_opts.lambda)*xA_mc
  z_mc = Rnd_Out_Z_All(z_mct,mc_opts.aff_correct_eps)

  p_mc::Vector{HybridMC{np,V,T}} = HybridMC{np,V,T}[HybridMC{np,V,T}(SMCg{np,V,T}(pmid[i],pmid[i],sone,sone,V(P[i].lo,P[i].hi),false)) for i=1:np]
  pref_mc::Vector{HybridMC{np,V,T}} = copy(p_mc)
  aff_mc::Vector{HybridMC{np,V,T}} = HybridMC{np,V,T}[HybridMC{np,V,T}(SMCg{np,V,T}(xA_mc[i].SMC.cc,xa_mc[i].SMC.cv,szero,szero,V(xa_mc[i].SMC.Intv.lo,xA_mc[i].SMC.Intv.hi),false)) for i=1:nxi]
  sto_out::Vector{Vector{HybridMC{np,V,T}}} = Vector{HybridMC{np,V,T}}[x_mc for j=1:(mc_opts.kmax+1)]
  sto_out[1] = copy(x_mc)
  optc = Any[szero,sone]

  for k=1:mc_opts.kmax
    PSMCg_Kernel!(h,hj,z_mc,aff_mc,p_mc,x_mc,mc_opts)
    Affine_Exp!(x_mc,p_mc,p_mc,xa_mc,xA_mc,z_mc,mc_opts)
    z_mc = Rnd_Out_Z_All(z_mc,mc_opts.aff_correct_eps)
    Correct_Exp!(z_mc,x_mc,X,nxi,np,mc_opts.aff_correct_eps)
    aff_mc = HybridMC{np,V,Float64}[HybridMC{np,V,Float64}(SMCg{np,V,Float64}(xA_mc[i].SMC.cc,
                                    xa_mc[i].SMC.cv,xA_mc[i].SMC.cc_grad,xa_mc[i].SMC.cv_grad,
                                    V(xA_mc[i].SMC.Intv.lo,xA_mc[i].SMC.Intv.hi),false)) for i=1:nxi]
    # store relaxation
    sto_out[k+1] = copy(x_mc)
    end
  return sto_out[end]
end

# Layout Concentration[1:9900], Pressure[9901:10000], Kt[10001], Lpt[10002]
function LowerBound(X::Vector{Interval{Float64}},k::Int64,pos::Int64,opt::Any,UBD::Float64,rs,cref,kd,D,N)

    rmax = 2
    conc_n = N*N

    opt[1].solver.SubGradRefine && set_hybrid_box!(X,x0,true)

    R = 1.0
    Svt = 200.0 # fix me
    Pvv = 1.0 # fix me
    Pv = 25.0
    lambda = 0.5
    opt = mc_opt(Float64)

    #Y = ?  # Setup initial relaxations
    sone::SVector{2,Float64} = ones(SVector{2,Float64})
    szero::SVector{2,Float64} = zeros(SVector{2,Float64})

    p_mc = HybridMC{2,Interval{Float64},Float64}[HybridMC{2,Interval{Float64},Float64}(
                                                 SMCg{2,Interval{Float64},Float64}(
                                                 mid.(X[i]),mid.(X[i]),sone,sone,Interval{Float64}(X[i].lo,X[i].hi),false)
                                                 ) for i=(N^2+N+1):(N^2+N+2)]
    pref = copy(p_mc)
    #Kt = p_mc[end-1]
    #Lpt = p_mc[end]

    # Get bounds for Kt and Lpt
    Perm,sigma = solutePerm(p_mc[end],rs)

    # Bound pressure via implicit form
    Pressure_M = Pressure_Matrix(p_mc[end-1],p_mc[end],N,R,Svt,Pvv)
    Pressure_RHS = Pressure_RHS(p_mc[end-1],p_mc[end],N,R,Svt,Pvv)
    midPress = [mid(Pressure_M[i,j].SMC.Intv) for i=1:N, j=1:N]
    FP = lufact(midPress)   # Gets LU factorization
    Pressure_M = FP\Pressure_M    # Precondition the matrix
    Pressure_RHS = FP\Pressure_RHS
    Press_mc::Vector{HybridMC{2,Interval{Float64},Float64}} = HybridMC{2,Interval{Float64},Float64}[HybridMC{2,Interval{Float64},Float64}(
                                                     SMCg{2,Interval{Float64},Float64}(
                                                     X[i].hi,X[i].lo,szero,szero,Interval{Float64}(X[i].lo,X[i].hi),false)
                                                     ) for i=(conc_n+1):(conc_n+N)]
    Press_bnds::Vector{Interval{Float64}} = [Press_mc.SMC.Intv for i=1:N]
    Press_xa_mc::Vector{HybridMC{2,Interval{Float64},Float64}} = HybridMC{2,Interval{Float64},Float64}[HybridMC{2,Interval{Float64},Float64}(
                                                                                                       SMCg{2,Interval{Float64},Float64}(X[i].lo,X[i].lo,szero,szero,Interval{Float64}(X[i].lo,X[i].lo),false)
                                                                                                       ) for i=(conc_n+1):(conc_n+N]
    Press_xA_mc::Vector{HybridMC{2,Interval{Float64},Float64}} = HybridMC{2,Interval{Float64},Float64}[HybridMC{2,Interval{Float64},Float64}(
                                                                                                       SMCg{2,Interval{Float64},Float64}(X[i].hi,X[i].hi,szero,szero,Interval{Float64}(X[i].hi,X[i].hi),false)
                                                                                                       ) for i=(conc_n+1):(conc_n+N]
    Press_z_mc::Vector{HybridMC{2,Interval{Float64},Float64}} = lambda*Press_xa_mc+(one(T)-lambda)*Press_xA_mc
    for r = 1:rmax
        SMCg_Dense_Newton_GS!(Press_z_mc,Press_mc,Pressure_M,Pressure_RHS,opt)     # Runs contractor                                                                   # Run Contractor
        Affine_Exp!(Press_mc,p_mc,p_mc,Press_xa_mc,Press_xA_mc,Press_z_mc,opt)     # Affine Expansion
        Correct_Exp!(Press_z_mc,Press_mc,Press_bnds,nx,np,0.0)      # Corrects if affine is out of bounds
    end

    # Calculate V bounds
    V = vcat(zero(T),diff(Press_mc))

    # Bound concentration via implicit form
    Conc_M = Concentration_Matrix(Press_mc,V,Peff,sigma,Svt,h,kd,D,N)
    midConc = [mid(Conc_M[i,j].SMC.Intv) for i=1:opt.nx, j=1:opt.nx]
    F = lufact(midConc)   # Gets LU factorization
    Conc_M = F\Conc_M    # Precondition the matrix

    #Conc_mc::Vector{HybridMC{2,Interval{Float64},Float64}} = ?  # set initial condition
    Conc_xa_mc::Vector{HybridMC{2,Interval{Float64},Float64}} = HybridMC{2,Interval{Float64},Float64}[HybridMC{2,Interval{Float64},Float64}(
                                                                                                       SMCg{2,Interval{Float64},Float64}(X[i].lo,X[i].lo,szero,szero,Interval{Float64}(X[i].lo,X[i].lo),false)
                                                                                                       ) for j=(N+1):(conc_n+N)]
    Conc_xA_mc::Vector{HybridMC{2,Interval{Float64},Float64}} = HybridMC{2,Interval{Float64},Float64}[HybridMC{2,Interval{Float64},Float64}(
                                                                                                       SMCg{2,Interval{Float64},Float64}(X[i].hi,X[i].hi,szero,szero,Interval{Float64}(X[i].hi,X[i].hi),false)
                                                                                                       ) for j=(1):(nxi)]
    Conc_z_mc::Vector{HybridMC{2,Interval{Float64},Float64}} = lambda*Conc_xa_mc+(one(Float64)-lambda)*Conc_xA_mc

    Conc_IntVal = [i<N ? 0.0 : 1.0  for i=1:N]
    Conc_bnds = Interval{Float64}[Interval(Conc_Intval[i] for i=1:N]
    Y = zeros(HybridMC{2,Interval{Float64},Float64},conc_n)
    Y[1:N] = Conc_IntVal
    Conc_mc = Y[1:N]
    for j = 1:N
        for r = 1:rmax
            #Conc_RHS = F\(Conc_RHS(i,Peff,sigma,Svt,h,kd,D,N) + Conc_mc   # Precondition righthand side
            #Conc_mc = HybridMC{2,Interval{Float64},Float64}[HybridMC{2,Interval{Float64},Float64}(SMCg{2,Interval{Float64},Float64}(X[i].lo,X[i].lo,szero,szero,Interval{Float64}(X[i].lo,X[i].lo),false)) for i=(N*(j-1)+1):(N*j)]
            #Conc_bnds = Interval{Float64}[Intv(Conc_mc) for i=1:N]
            #Conc_xa_mc = HybridMC{2,Interval{Float64},Float64}[HybridMC{2,Interval{Float64},Float64}(SMCg{2,Interval{Float64},Float64}(X[i].lo,X[i].lo,szero,szero,Interval{Float64}(X[i].lo,X[i].lo),false)) for i=(N*(j-1)+1):(N*j)]
            #Conc_xA_mc = HybridMC{2,Interval{Float64},Float64}[HybridMC{2,Interval{Float64},Float64}(SMCg{2,Interval{Float64},Float64}(X[i].hi,X[i].hi,szero,szero,Interval{Float64}(X[i].hi,X[i].hi),false)) for i=(N*(j-1)+1):(N*j)]
            #Conc_z_mc = lambda*Conc_xa_mc + (1.0-lambda)*Conc_xA_mc
            #SMCg_Dense_Newton_GS!(Conc_z_mc,Conc_mc,Conc_M,Conc_RHS,opt)                                                      # Runs contractor                                                                   # Run Contractor
            #Affine_Exp!(Conc_mc,p_mc,p_mc,Conc_xa_mc,Conc_xA_mc,Conc_z_mc,opt)                                               # Affine Expansion
            #Correct_Exp!(Conc_z_mc,Conc_mc,Conc_bnds,nx,np,0.0)                                             # Corrects if affine is out of bounds
            #Y[(N*(j-1)+1):(N*j)] = Conc_mc
            # Stores current concentraton values
        end
    end

    # calculate the SSE for the objective
    obj = 0.0
    for i=2:N
        obj = obj  + (cref[i-1]-mean(Y[(N*(i-1)+1):(N*(i-1)+N)]))^2
    end

    # return the object function
    return obj
end

function UpperBound((X::Vector{Interval{Float64}},k::Int,pos::Int,opt,temp,rs,cref,kd,D,N)

    # Unpack variables
    R = 1.0
    Svt = 200.0
    Pvv = 1.0
    N = 100
    Y = mid.(X[1:(end-2)])
    Lpt = mid.(X[end-1])

    # Deen Calculations
    Perm,sigma = solutePerm(mid.(Lpt,rs)
    Kt = mid.(X[end])

    # calculate the pressure and velocity field
    Pressure_M = Pressure_Matrix(Lpt,Kt,N,R,Svt,Pvv)
    Pressure_RHS = Pressure_RHS(Lpt,Kt,N,R,Svt,Pvv)
    P = Pressure_M\Pressure_RHS
    V = vcat(zero(T),diff(P))

    # calculate the concentration field
    Conc_M = Concentration_Matrix(P,V,Peff,sigma,Svt,h,kd,D,N)
    inv_M = inv(Conc_M)
    M = inv_M*Conc_M
    rhs = Conc_RHS(i+1)
    for i=1:(N-1)
        Conc_RHS = inv_M*(rhs + mid.(Y[(N*(i-1)+1):(N*(i-1)+N)]))
        Y[(N*i+1):(N*i+N)] = M\Conc_RHS
    end

    # calculate the SSE for the objective
    obj = 0.0
    for i=2:N
        obj = (cref[i-1]-mean(Y[(N*(i-1)+1):(N*(i-1)+N)]))^2
    end

    # return the object function
    return obj
end
