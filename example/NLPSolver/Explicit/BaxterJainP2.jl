workspace()
using EAGO
using JuMP
using Ipopt
using MathProgBase
using BenchmarkTools
using IntervalArithmetic
using Plots

function Isolated_Pressure(N::Int,R::Float64,Lpt::T,Svt::Float64,Kt::Float64,Pvv::Float64) where T<:Real

    # Pressure Profile
    att = R*sqrt(Lpt*Svt/Kt)

    # N -- number of spatial grid points
    r = linspace(0,R,N)
    r = r/R
    dr = 1/(N-1)
    ivdr2 = 1./dr^2
    att2 = att^2
    A = zeros(N,N)
    F = zeros(N,1)

    #at = zeros(N,1)

    M = N

    #for i in 1:M
#        at[i] = att
#    end

    for i in 2:M-1
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
    F[1] = -att2*Pvv
    A[N,N] = 1.0

    # Solution of linear problem for pressure distribution
    P = A\F
    return P
end

function Isolated_Pressure_Form(N::Int,R::Float64,Lpt::T,Svt::Float64,Kt::T,Pvv::Float64) where T<:Real

    att = R*sqrt(Lpt*Svt/Kt)
    #println("att: $att")
    r = Vector(linspace(0,R,N))/R
    dimP = zeros(T,N)
    dimV = zeros(T,N)
    dimP[2:end] = one(T) - sinh.(att*r[2:end])./(sinh(att)*r[2:end])
    #dimP[2:end] = [0.01*one(T) for i=2:N]
    dimV[2:end] = (att*r[2:end].*cosh.(att*r[2:end])-sinh.(att*r[2:end]))./(sinh(att)*r[2:end].^2)
    #dimV[2:end] = [one(T) for i=2:N]

    Pinf = zero(T)
    Pdel = Pvv - Pinf
    P = Pdel*dimP+Pinf
    V = (Kt*Pdel/R)*dimV
    P[1] = one(T)

    return P,V
end

function MST_Form(t::Float64,c::Vector{T},P::Vector{T},V::Vector{T},N::Int,sigma::T,
                Peff::T,Lpt::T,Svt::Float64,Kt::T,Pvv::Float64,
                Pv::Float64,D::Float64,r::Vector{Float64},dr::Float64,kd::Float64) where T<:Real

    co::Float64 = 1.0 # dimensionless drug concentration
    tspan::Float64 = 1.0*3600.0 # length of simulation type

    f::Vector{T} = zeros(T,N)
    cv::Float64 = co*exp(-t/kd/3600.0)  # vascular concentration of the drug following exponential decay
    coeff1::T = Peff*Svt
    coeff2::T = Lpt*(Svt*Pv*cv*(1.0-sigma))
    coeff3::Float64 = 2.0*D/dr
    coeff4::Float64 = D/dr^2
    f[1] = 2.0*D*(c[2]-c[1])/dr^2 + Peff*Svt*(cv-c[1]) + coeff2*(Pvv-P[1])

    for j in 2:N-1
        f[j] = ((coeff3/r[j])*((c[j+1]-c[j])) + coeff4*(c[j+1]-2.0*c[j]+c[j-1]) +
               V[j]*((c[j+1]-c[j])/dr) + coeff1*(cv-c[j])) + coeff2*(Pvv-P[j])
    end

    f[N] = zero(T)
    return f
end

function RK4_Form(x0::Float64,y0::Vector{T},h::Float64,n_out::Int,i_out::Int,
                  P::VecOrMat{T},V::Vector{T},N::Int64,sigma::T,Peff::T,Lpt::T,Svt::Float64,
                  Kt::T,Pvv::Float64,Pv::Float64,D::Float64,r::Vector{Float64},
                  dr::Float64,kd::Float64) where T<:Real

    xout = zeros(T,n_out+1)
    yout = zeros(T,n_out+1,length(y0))
    xout[1] = x0
    yout[1,:] = y0
    x = x0
    y = y0
    for j = 2:n_out+1
        for k = 1:i_out
            k1 = MST_Form(x,y,P,V,N,sigma,Peff,Lpt,Svt,Kt,Pvv,Pv,D,r,dr,kd)
            k2 = MST_Form(x+0.5*h,y+0.5*h*k1,P,V,N,sigma,Peff,Lpt,Svt,Kt,Pvv,Pv,D,r,dr,kd)
            k3 = MST_Form(x+0.5*h,y+0.5*h*k2,P,V,N,sigma,Peff,Lpt,Svt,Kt,Pvv,Pv,D,r,dr,kd)
            k4 = MST_Form(x,y+h*k3,P,V,N,sigma,Peff,Lpt,Svt,Kt,Pvv,Pv,D,r,dr,kd)
            y = min(max(y + (h/6.0)*(k1 + 2.0*k2 + 2.0*k3 + k4),0.0),1.0)
            #y = y + (h/6.0)*(k1 + 2.0*k2 + 2.0*k3 + k4)
            x = x + h
        end
        xout[j] = x
        yout[j,:] = y
    end
    return xout, yout
end

function EE_Form(x0::Float64,y0::Vector{T},h::Float64,n_out::Int,i_out::Int,
                  P::VecOrMat{T},V::Vector{T},N::Int64,sigma::Float64,Peff::Float64,Lpt::T,Svt::Float64,
                  Kt::Float64,Pvv::Float64,Pv::Float64,D::Float64,r::Vector{Float64},
                  dr::Float64,kd::Float64) where T<:Real

      xout = zeros(T,n_out+1)
      yout = zeros(T,n_out+1,length(y0))
      xout[1] = x0
      yout[1,:] = y0
      x = x0
      y = y0
      for j = 2:n_out+1
          for k = 1:i_out
              y = max(0.0,min(y + h*MST_Form(x,y,P,V,N,sigma,Peff,Lpt,Svt,Kt,Pvv,Pv,D,r,dr,kd),1.0))
              x = x + h
          end
          xout[j] = x
          yout[j,:] = y
      end
      return xout, yout
end


function Isolated_Model_Form(N::Int,Kt::Float64,Lpt::T,Svt::Float64,D::Float64,
                        sigma::Float64,Peff::Float64,R::Float64,Pv::Float64,
                        Pvv::Float64,kd::Float64,n_nodes::Int) where T<:Real

    r = Vector(linspace(0,R,N))
    r = r/R
    dr = 1.0/(N-1)

    # Solution of steady state pressure model
    P,V = Isolated_Pressure_Form(N,R,Lpt,Svt,Kt,Pvv)

    # Initial solute concentration
    c_0 = zeros(T,N)
    #c_0[N] = one(T)

    time_end = 1*3600.0 # length of simulation (seconds)
    n_out = n_nodes - 1
    h = time_end/n_out;
    i_out = 1

    time, c = RK4_Form(0.0,c_0,h,n_out,i_out,P,V,N,sigma,Peff,Lpt,Svt,Kt,Pvv,Pv,D,r,dr,kd)
    return time, c
end

function Fitting_Objective(N::Int,Kt::T,Lpt::T,Svt::Float64,D::Float64,
                           sigma::T,Peff::T,R::Float64,Pv::Float64,
                           Pvv::Float64,kd::Float64,n_nodes::Int,n_time::Int,
                           cref::VecOrMat{Float64}) where T<:Real

    r = Vector(linspace(0,R,N))
    r = r/R
    dr = 1.0/(N-1)

    # Solution of steady state pressure model
    P,V = Isolated_Pressure_Form(N,R,Lpt,Svt,Kt,Pvv)

    # Initial solute concentration
    c_0 = zeros(T,N)
    #c_0[N] = one(T)

    time_end = 1*3600.0 # length of simulation (seconds)
    n_out = n_nodes - 1
    h = time_end/n_out;
    i_out = 1

    time, c = RK4_Form(0.0,c_0,h,n_out,i_out,P,V,N,sigma,Peff,Lpt,Svt,Kt,Pvv,Pv,D,r,dr,kd)
    #println("c[100,30]: $(c[100,30])")
    #println("c: $c")
    SSE = zero(T)
    for j = 1:(n_time)
        c_model = mean(c[j*Int(floor(n_nodes/n_time)),:])
        #println("c_model vec at $j: $(c[j*Int(floor(n_nodes/n_time)),:])")
        #println("c_model at $j: $c_model")
        SSE = SSE + (c_model-cref[j])^2
    end

    return SSE
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
    D0 = kB*T/(6*pi*etac*r_partc)*1.0e4;   # Diffusivity (cm^2/s)

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

function Fit_With_Deen(N,x,Svt,D,R,Pv,Pvv,kd,rs,n_nodes,n_time,cref)
    Perm,sigma = solutePerm(x[1],rs)
    output = Fitting_Objective(N,x[2],x[1],Svt,D,sigma,Perm,R,Pv,Pvv,kd,n_nodes,n_time,cref)
    return output
end

idx = 2

# Input space nodes and time nodes
#N = 20
#n_nodes = 20 #6000
#n_time = 20

N = 50
n_nodes = 50 #6000
n_time = 50

# Parameters for creating data for fit
co = 1
t = (3600/n_time).*Array(1:n_time)
s = EAGO_NLPSolver(LBD_func_relax = "NS-STD-OFF",
                   LBDsolvertype = "LP",
                   LBD_func_relax = "Interval",
                   LBDsolvertype = "Interval",
                   UBDsolvertype = "Interval",
                   #UBDsolvertype = "Ipopt",
                   probe_depth = -1,
                   variable_depth = -1,
                   DAG_depth = -1,
                   STD_RR_depth = -1,
                   validated = true)

# Input known metabolic parameters for Rhodamine
Peff_set = [9.59873e-07;4.6098e-06;2.80047e-06]
#Kt = 0.9e-7              # Hydraulic conductivity of tumor
#Lpt = 5.0e-7             # hydraulic conductivity of tumor vessels
Svt = 200.0                # tumor vascular density
D = 2.0e-7                 # solute diffusion coefficient
#sigma = 0.0                # solute reflection coefficient
#Perm = Peff_set[idx]*0.8 # solute vascular permeability
R = 1.0                   # tumor radius (cm)
Pv = 25.0                  # vascular pressure (mmHg)
Pvv = 1.0                 # vascular pressure dimensionless
kd = 1480*60.0             #blood circulation time of drug in hours
rs = 13/2.0
#=
a = zeros(100)
for i=1:100
    j = (i-1)*((6e-7)-(6e-9))/(100.0)+6e-9
    a[i] = Fitting_Objective(N,Kt,j,Svt,D,sigma,Peff_set[1],R,Pv,Pvv,kd,n_nodes,n_time,cref1)
end
=#
#plotly()
#plot(a)
#gui()

#@btime Fitting_Objective(N,Kt,6e-9,Svt,D,sigma,Peff_set[1],R,Pv,Pvv,kd,n_nodes,n_time,cref1)
#@btime Fitting_Objective(N,Kt,6e-7,Svt,D,sigma,Peff_set[1],R,Pv,Pvv,kd,n_nodes,n_time,cref1)
#plotly()
#plot(tcalc,ccalc[:,10])
#plot!(tcalc1,ccalc1[:,10])
#gui()
#@time Isolated_Model_Form(N,Kt,Lpt,Svt,D,sigma,Peff_set[1],R,Pv,Pvv,kd,n_nodes)
#IIval = Fitting_Objective(N,Kt,Lpt,Svt,D,sigma,Peff_set[1],R,Pv,Pvv,kd,n_nodes,n_time,cref1)
#@time Fitting_Objective(N,Kt,Interval(6e-9,6e-7),Svt,D,sigma,Peff_set[1],R,Pv,Pvv,kd,n_nodes,n_time,cref1)

cref1 = (co*Peff_set[1]*Svt*kd/(1-Peff_set[1]*Svt*kd))*(exp(-Peff_set[1]*Svt*t)-exp(-t/kd))
cref2 = (co*Peff_set[2]*Svt*kd/(1-Peff_set[2]*Svt*kd))*(exp(-Peff_set[2]*Svt*t)-exp(-t/kd))
cref3 = (co*Peff_set[3]*Svt*kd/(1-Peff_set[3]*Svt*kd))*(exp(-Peff_set[3]*Svt*t)-exp(-t/kd))

#=
println("IntervalEval 1")
IntvObj = Fitting_Objective(N,Kt,Interval(6e-19,6e-1),Svt,D,sigma,Peff_set[1],R,Pv,Pvv,kd,n_nodes,n_time,cref1)
println("IntervalEval 2")
IntvObj1 = Fitting_Objective(N,Kt,Interval(6e-9,6e-7),Svt,D,sigma,Peff_set[2],R,Pv,Pvv,kd,n_nodes,n_time,cref2)
println("IntervalEval 3")
IntvObj2 = Fitting_Objective(N,Kt,Interval(6e-9),Svt,D,sigma,Peff_set[2],R,Pv,Pvv,kd,n_nodes,n_time,cref2)
IntvObj3 = Fitting_Objective(N,Kt,Interval(6e-9,1e-8),Svt,D,sigma,Peff_set[1],R,Pv,Pvv,kd,n_nodes,n_time,cref1)
IntvObj4 = Fitting_Objective(N,Kt,Interval(6e-9,6.01e-9),Svt,D,sigma,Peff_set[1],R,Pv,Pvv,kd,n_nodes,n_time,cref1)
IntvObj5 = Fitting_Objective(N,Kt,Interval(6e-9,6.01e-9),Svt,D,sigma,Peff_set[1],R,Pv,Pvv,kd,n_nodes,n_time,cref1)
IntvObj6 = Fitting_Objective(N,Kt,Interval(6e-9,6.0001e-9),Svt,D,sigma,Peff_set[1],R,Pv,Pvv,kd,n_nodes,n_time,cref1)
IntvObj7 = Fitting_Objective(N,Kt,Interval(6e-9,6.000001e-9),Svt,D,sigma,Peff_set[1],R,Pv,Pvv,kd,n_nodes,n_time,cref1)
IntvObj8 = Fitting_Objective(N,Kt,Interval(6e-9,6.00000001e-9),Svt,D,sigma,Peff_set[1],R,Pv,Pvv,kd,n_nodes,n_time,cref1)
IntvObj9 = Fitting_Objective(N,Kt,Interval(6.000005e-9,6.0000051e-9),Svt,D,sigma,Peff_set[1],R,Pv,Pvv,kd,n_nodes,n_time,cref1)
#IntvObj = Fitting_Objective(N,Kt,Interval(6e-9,6e-7),Svt,D,sigma,Peff_set[1],R,Pv,Pvv,kd,n_nodes,n_time,cref1)
IntvP,IntvV = Isolated_Pressure_Form(N,R,6e-9,Svt,Kt,Pvv)
=#

# Fits the data for control, rhodamine
#tcalc, ccalc = Isolated_Model_Form(N,Kt,Lpt,Svt,D,sigma,Peff_set[1],R,Pv,Pvv,kd,n_nodes)
#objv = Fitting_Objective(N,Kt,Lpt,Svt,D,sigma,Peff_set[1],R,Pv,Pvv,kd,n_nodes,n_time,cref1)

function local_solve(s,X)
    cs = EAGO.callback_storage()
    cs.IPOPT_UBD_eval_grad_f! = (x::Vector{Float64}, f_grad::Vector{Float64}) -> EAGO.IPOPT_UBD_eval_grad_f!(x, f_grad, s.Opts)
    cs.IPOPT_UBD_eval_g! = (x::Vector{Float64}, g::Vector{Float64}) -> EAGO.IPOPT_UBD_eval_g!(x, g, s.Opts)
    cs.IPOPT_UBD_eval_jac_g! = (x::Vector{Float64}, mode::Symbol, rows::Vector{Int32}, cols::Vector{Int32}, values::Array{Float64,1}) -> EAGO.IPOPT_UBD_eval_jac_g!(x, mode, rows, cols, values, s.Opts, cs)
    cs.IPOPT_UBD_eval_h = (x::Vector{Float64}, mode::Symbol,rows::Vector{Int32}, cols::Vector{Int32}, obj_factor::Float64, lambda::Vector{Float64}, values::Array{Float64,1}) -> EAGO.IPOPT_UBD_eval_h(x, mode, rows, cols, obj_factor, lambda, values, s.Opts)
    cols_temp = zeros(Int32,s.Opts.numVar*s.Opts.numConstr)
    rows_temp = zeros(Int32,s.Opts.numVar*s.Opts.numConstr)
    for i = 1:s.Opts.numConstr
        cols_temp[(s.Opts.numVar*(i-1)+1):(s.Opts.numVar*i)] = 1:s.Opts.numVar
        rows_temp[(s.Opts.numVar*(i-1)+1):(s.Opts.numVar*i)] = ones(s.Opts.numVar)*i
    end
    cs.col_temp_Ipopt_LBD = cols_temp
    cs.row_temp_Ipopt_LBD = rows_temp

    val, pnt, feas, outer = EAGO.Ipopt_UBD(X,1,1,[s.Opts,cs],[])
    return val, pnt, feas
end

# partition each dimension into r parts and solve (grid multi-start)
function multi_start(s,X,r)
    vmin = Inf
    pmin = mid.(X)
    for i=1:r^2
        #Xt = [X[]]
        v,p,f = local_solve(s,X)
        if ((v < vmin) && f)
            vmin = v
            pmin = p
        end
    end
    return vmin,pmin
end

f1(x) = Fit_With_Deen(N,x,Svt,D,R,Pv,Pvv,kd,rs,n_nodes,n_time,cref1)
m1 = MathProgBase.NonlinearModel(s)
MathProgBase.loadproblem!(m1, 2, 0, [1.0e-7,1.0e-7], [1.2e-6,0.8e-6],[], [], :Min, f1, [])
X = [Interval(1.0e-7,1.2e-6), Interval(1.0e-7,0.8e-6)]
ls = local_solve(m1,X)
ms = multi_start(m1,X,1)
#=
MathProgBase.optimize!(m1)

# Fits the data for 3mg/kg, rhodamine
f2(x) = Fit_With_Deen(N,x,Svt,D,R,Pv,Pvv,kd,rs,n_nodes,n_time,cref2)
m2 = MathProgBase.NonlinearModel(s)
MathProgBase.loadproblem!(m2, 2, 0, [1.0e-7,0.5e-6], [1.5e-6,2.5e-6],[], [], :Min, f2, [])
MathProgBase.optimize!(m2)

# Fits the data for 30mg/kg, rhodamine
f3(x) = Fit_With_Deen(N,x,Svt,D,R,Pv,Pvv,kd,rs,n_nodes,n_time,cref3)
m3 = MathProgBase.NonlinearModel(s)
MathProgBase.loadproblem!(m3, 2, 0, [1.0e-7,1.4e-6], [1.0e-6,2.0e-6],[], [], :Min, f3, [])
MathProgBase.optimize!(m3)


# Input known metabolic parameters
N = 50
n_nodes = 50 #6000
n_time = 50

Peff_set = [8.18378e-07;4.30307e-06;1.62231e-06]
#Kt = 0.9e-7              # Hydraulic conductivity of tumor
#Lpt = 5.0e-6             # hydraulic conductivity of tumor vessels
Svt = 200.0                # tumor vascular density
D = 1.4375e-07                 # solute diffusion coefficient
#sigma = 0.0                # solute reflection coefficient
#Perm = Peff_set[idx]*0.8 # solute vascular permeability
R = 1.0                   # tumor radius (cm)
Pv = 25.0                  # vascular pressure (mmHg)
Pvv = 1.0                 # vascular pressure dimensionless
kd = 1278*60.0             #blood circulation time of drug in hours
rs = 32/2

# Fits the data for control, rhodamine
cref4 = (co*Peff_set[1]*Svt*kd/(1-Peff_set[1]*Svt*kd))*(exp(-Peff_set[1]*Svt*t)-exp(-t/kd))
f4(x) = Fit_With_Deen(N,x,Svt,D,R,Pv,Pvv,kd,rs,n_nodes,n_time,cref4)
m4 = MathProgBase.NonlinearModel(s)
MathProgBase.loadproblem!(m4, 2, 0, [2.0e-8,1.0e-6], [1.0e-7,0.8e-6],[], [], :Min, f4, [])
MathProgBase.optimize!(m4)

# Fits the data for 3mg/kg, rhodamine
cref5 = (co*Peff_set[2]*Svt*kd/(1-Peff_set[2]*Svt*kd))*(exp(-Peff_set[2]*Svt*t)-exp(-t/kd))
f5(x) = Fit_With_Deen(N,x,Svt,D,R,Pv,Pvv,kd,rs,n_nodes,n_time,cref5)
m5 = MathProgBase.NonlinearModel(s)
MathProgBase.loadproblem!(m5, 2, 0, [1.0e-7,1.5e-6], [1.5e-6,2.5e-6],[], [], :Min, f5, [])
MathProgBase.optimize!(m5)

# Fits the data for 30mg/kg, rhodamine
cref6 = (co*Peff_set[3]*Svt*kd/(1-Peff_set[3]*Svt*kd))*(exp(-Peff_set[3]*Svt*t)-exp(-t/kd))
f6(x) = Fit_With_Deen(N,x,Svt,D,R,Pv,Pvv,kd,rs,n_nodes,n_time,,cref6)
m6 = MathProgBase.NonlinearModel(s)
MathProgBase.loadproblem!(m6, 2, 0, [1.0e-7,1.4e-6], [1.5e-6,2.5e-6],[], [], :Min, f6, [])
MathProgBase.optimize!(m6)


tcalc, ccalc = Isolated_Model_Form(N,Kt,Lpt,Svt,D,sigma,Peff_set[1],R,Pv,Pvv,kd,n_nodes)
tcalc1, ccalc1 = Isolated_Model_Form(N,Kt,Lpt,Svt,D,sigma,Peff_set[1],R,Pv,Pvv,kd,n_nodes)
objv1 = Fitting_Objective(N,Kt,6e-9,Svt,D,sigma,Peff_set[1],R,Pv,Pvv,kd,n_nodes,n_time,cref1)
objv2 = Fitting_Objective(N,Kt,6e-7,Svt,D,sigma,Peff_set[1],R,Pv,Pvv,kd,n_nodes,n_time,cref1)

=#
