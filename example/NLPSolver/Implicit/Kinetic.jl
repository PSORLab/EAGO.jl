
###############################################################################
# Test Problem #0 (Kinetic Data Fitting Problem)
###############################################################################
# LOADS DATA
I_Data = [66.0952, 104.762, 110.333, 114.905, 122.238, 125.429, 125.429, 123.476, 121.286,
          118.857, 117.667, 116.143, 113.857, 111.571, 108.81,  105.952, 104.048, 102.048,
          100.143, 98.5238, 96.2381, 94.381,  91.6667, 89.5714, 87.1429, 84.8571, 83.4286,
          81.1905, 78.9048, 77.0476, 75.4762, 73.4762, 71.8095, 70.6667, 68.381,  67.3333,
          65.0952, 63.7143, 62.0476, 60.8571, 59.619,  58.2857, 57.4762, 56.4762, 55.8095,
          54.5238, 53,      51.8571, 50.4286, 49.381,  47.9524, 47.3714, 46.8952, 46.4857,
          45.9048, 45.0762, 44.3238, 43.4143, 43.5429, 42.3619, 41.8381, 40.2381, 39.1286,
          38.7857, 37.081,  36.9524, 36.581,  36.281,  35.3476, 34.8905, 34.1667, 33.6714,
          32.9667, 31.8429, 31.5429, 31.1476, 30.9905, 29.9571, 29.1333, 28.7857, 28.4429,
          28.3476, 27.5429, 27.4333, 27.6048, 27.1762, 27.2,    26.4333, 25.7619, 24.8095,
          24.7429, 24.2857, 24.1714, 23.5667, 23.5476, 23.3952, 22.919,  22.3095, 21.8048,
          21.2857, 21.2048, 20.8429, 20.4429, 20.0048, 19.9381, 19.5,    19.8667, 18.9333,
          19.1381, 18.9619, 18.5476, 17.9048, 17.7571, 18.5333, 18.3762, 18.3571, 18.3286,
          18.2762, 18.3952, 17.5952, 18.1524, 18.1952, 17.8476, 17.9095, 17.5048, 17.5,
          15.9619, 16.2095, 16.181,  15.6952, 15.7095, 15.4619, 15.9476, 16,      16.1952,
          16.1143, 15.7429, 15.5762, 15.7048, 15.8095, 15.6667, 14.9048, 14.5857, 14.7524,
          14.7571, 14.9762, 14.5333, 14.5524, 14.0143, 13.6286, 13.4429, 13.4667, 13.319,
          12.9333, 13.1238, 12.7476, 12.9333, 13.0714, 13.0714, 12.7619, 12.4238, 12.5143,
          12.9143, 12.5714, 13.3667, 13.2286, 13.7905, 13.7571, 13.5905, 12.9667, 12.981,
          12.8857, 12.919,  13.0143, 13.0095, 12.3857, 12.5571, 12.3429, 12.7571, 12.681,
          12.5429, 12.1857, 12.7905, 12.5571, 12.8429, 12.5476, 12.5714, 12.3762, 11.9952,
          11.4571, 11.3,    11.1524, 11.681,  11.619,  11.9048, 12,      12.0762, 11.9143,
          11.7619, 11.5333]
n = 60
#n = 200

# SETS UP OBJECTIVE FUNCTION FOR IMPLICIT SOLVER
function fi(z,p,I)
    n = 60
    #n = 200
    Icalc =  zeros(typeof(p[1]),n)
    for i=1:n
        Icalc[i] = z[5*(i-1)+1] + (2.0/21.0)*z[5*(i-1)+2] + (2.0*21.0)*z[5*(i-1)+3]
    end
    return sum((Icalc-I).^2)
end
f(z,p) = fi(z,p,I_Data[1:n])

# SETS UP CONSTRAINT FUNCTION FOR IMPLICIT SOLVER
function h(out,z,p)
    #println("started h eval")
    n = 60 # n = 200
    delT = 0.01
    T = 273.0
    K2 = 46.0*exp(6500.0/T-18.0)
    K3 = 2*K2
    k1 = 53.0
    ks1 = k1*1.0E-6
    k5 = 1.2E-3
    cO2 = 2.0E-3
    k2f = p[1]
    k3f = p[2]
    k4 = p[3]

    #out = Array{typeof(p[1])}(5*n)
    out[1] = -z[1] + delT*(k1*z[4]*z[5] - cO2*(k2f+k3f)*z[1] + (k2f/K2)*z[3] + (k3f/K3)*z[2] - k5*z[1]^2)
    out[2] = -z[2] + delT*(k3f*cO2*z[1]-(k3f/K3+k4)*z[2])
    out[3] = -z[3] + delT*(k2f*cO2*z[1]-(k2f/K2)*z[3])
    out[4] = 0.4 - z[4] + delT*(-k1*z[4]*z[5])
    out[5] = 140.0 - z[5] + delT*(-k1*z[4]*z[5])
    for i=2:n
         out[5*(i-1)+1] = z[5*(i-1)-4] - z[5*(i-1)+1] + delT*(k1*z[5*(i-1)+4]*z[5*(i-1)+5] - cO2*(k2f+k3f)*z[5*(i-1)+1] +
                                                             (k2f/K2)*z[5*(i-1)+3] + (k3f/K3)*z[5*(i-1)+2] - k5*z[5*(i-1)+1]^2)
         out[5*(i-1)+2] = z[5*(i-1)-3] - z[5*(i-1)+2] + delT*(k3f*cO2*z[5*(i-1)+1]-(k3f/K3+k4)*z[5*(i-1)+2])
         out[5*(i-1)+3] = z[5*(i-1)-2] - z[5*(i-1)+3] + delT*(k2f*cO2*z[5*(i-1)+1]-(k2f/K2)*z[5*(i-1)+3])
         out[5*(i-1)+4] = z[5*(i-1)-1] - z[5*(i-1)+4] + delT*(-k1*z[5*(i-1)+4]*z[5*(i-1)+5])
         out[5*(i-1)+5] = z[5*(i-1)] - z[5*(i-1)+5] + delT*(-k1*z[5*(i-1)+4]*z[5*(i-1)+5])
    end
    #println("end h eval")
    #return out
end

function hj(out,z,p)
    #println("started hj eval")
    n = 60 # n = 200
    delT = 0.01
    T = 273.0
    K2 = 46.0*exp(6500.0/T-18.0)
    K3 = 2*K2
    k1 = 53.0
    ks1 = k1*1.0E-6
    k5 = 1.2E-3
    cO2 = 2.0E-3
    k2f = p[1]
    k3f = p[2]
    k4 = p[3]

    #out = zeros(typeof(p[1]),5*n,5*n)
    out[1,1] = -1.0 - delT*(cO2*(k2f+k3f) + 2*k5*z[1])
    out[1,2] = delT*(k3f/K3)
    out[1,3] = delT*(k2f/K2)
    out[1,4] = delT*k1*z[5]
    out[1,5] = delT*k1*z[4]

    out[2,1] = delT*(k3f*cO2)
    out[2,2] = -1.0 - delT*(k3f/K3+k4)

    out[3,1] = delT*(k2f*cO2)
    out[3,3] = -1.0 - delT*(k2f/K2)

    out[4,4] = -1.0 + delT*(-k1*z[5])
    out[4,5] = delT*(-k1*z[4])

    out[5,4] = delT*(-k1*z[5])
    out[5,5] = -1.0 + delT*(-k1*z[4])

    for i=2:n
         out[5*(i-1)+1,5*(i-1)-4] = 1.0
         out[5*(i-1)+1,5*(i-1)+1] = -1.0 + delT*(-cO2*(k2f+k3f) - 2*k5*z[5*(i-1)+1])
         out[5*(i-1)+1,5*(i-1)+2] = delT*(k3f/K3)
         out[5*(i-1)+1,5*(i-1)+3] = delT*(k2f/K2)
         out[5*(i-1)+1,5*(i-1)+4] = delT*(k1*z[5*(i-1)+5])
         out[5*(i-1)+1,5*(i-1)+5] = delT*(k1*z[5*(i-1)+4])

         out[5*(i-1)+2,5*(i-1)-3] = 1.0
         out[5*(i-1)+2,5*(i-1)+2] = -1.0 + delT*(-(k3f/K3+k4))
         out[5*(i-1)+2,5*(i-1)+1] = delT*(k3f*cO2)

         out[5*(i-1)+3,5*(i-1)-2] = 1.0
         out[5*(i-1)+3,5*(i-1)+1] = delT*(k2f*cO2)
         out[5*(i-1)+3,5*(i-1)+3] = -1.0 + delT*(-(k2f/K2))

         out[5*(i-1)+4,5*(i-1)-1] = 1.0
         out[5*(i-1)+4,5*(i-1)+4] = -1.0 + delT*(-k1*z[5*(i-1)+5])
         out[5*(i-1)+4,5*(i-1)+5] = delT*(-k1*z[5*(i-1)+4])

         out[5*(i-1)+5,5*(i-1)] = 1.0
         out[5*(i-1)+5,5*(i-1)+5] = -1.0 + delT*(-k1*z[5*(i-1)+4])
         out[5*(i-1)+5,5*(i-1)+4] = delT*(-k1*z[5*(i-1)+5])

    end
    #println("ended hj eval")
    #return out
end

# SETS UP FUNCTIONS FOR EXPLICIT UPPER PROBLEM
fx = x -> f(x[1:(5*n)],x[((5*n)+1):((5*n)+3)])
#hx = x -> h(x[1:(5*n)],x[((5*n)+1):((5*n)+3)])
function hxi(x)
       n = 60
       temp = copy(x[1:5n])
       h(temp,x[1:(5*n)],x[((5*n)+1):((5*n)+3)])
       return temp
end
# SETS UP VARIABLE BOUNDS
YL = zeros(Float64,((5*n)+3))
YL[((5*n)+1)] = 10.0
YL[((5*n)+2)] = 10.0
YL[((5*n)+3)] = 0.001

YU = 140.0*ones(Float64,((5*n)+3))
YU[4:5:(5*n)] = 0.4
YU[((5*n)+1)] = 1200.0
YU[((5*n)+2)] = 1200.0
YU[((5*n)+3)] = 40.0

println("fx start")
f1 = fx(YL)
println("fx checks out")
#f2 = hx(YL)
println("hx checks out")
#f3 = hj(ones(5*n),ones(3))
println("hj checks out")
#forwardJac = (x,p) -> ForwardDiff.jacobian(hx,vcat(x,p))[1:5n,1:5n]

# SETS UP VALUE OF DOUBLESIDED CONSTRAINTS
hzero = zeros(5*n)

# SETS UP IMPLICIT OPTIONS FOR SOLVER
opt1 = Any[1    #Int64: Number of iterations
       1.0E-6 #Float64: Tolerance for equality of
       1.0E-6 #Float64: Add Interval(1E-8,1E8) to add to M[i,i] when
              #         processing extended interval division.
      ]

imp_opt = ImplicitSolver()
imp_opt.opts.style = "KrawczykCW" #"KrawczykCW" # sets MC contractor style
imp_opt.h = h
imp_opt.hj = hj
imp_opt.f = f
imp_opt.g = [] # inherites bounds from jump model
imp_opt.nx = 5*n
imp_opt.flag = true
imp_opt.Intv_Cntr =  "KrawczykCW" #"NewtonGS" # sets interval style
imp_opt.ParamInt = opt1
imp_opt.numConstr = 0
imp_opt.gL_Loc = []
imp_opt.gU_Loc = []
imp_opt.gL = []
imp_opt.gU = []
imp_opt.Inplace = true


function LBDKinetic(i,n)
  if (i == (5*n+1))
    return 10.0
elseif (i == (5*n+2))
    return 10.0
elseif (i == (5*n+3))
    return 0.001
  else
    return 0.0
  end
end

function UBDKinetic(i,n)
  if (i == (5*n+1))
    return 1200.0
elseif (i == (5*n+2))
    return 1200.0
elseif (i == (5*n+3))
    return 40.0
elseif (mod(i,5) == 4)
    return 0.4
  else
    return 140.0
  end
end

LBDKinetic(1,n)
UBDKinetic(2,n)

println("began param test")
Xtest1 = [Interval(LBDKinetic(i,n),UBDKinetic(i,n)) for i=1:(5*n)]
println("began param test 1")
Ptest1 = [Interval(LBDKinetic(i,n),UBDKinetic(i,n)) for i=(5*n+1):(5*n+3)]
println("began param test 2")
pmidtest1 = mid.([Interval(LBDKinetic(i,n),UBDKinetic(i,n)) for i=(5*n+1):(5*n+3)])
println("began param test 3")
println("imp_opt.opts: $(imp_opt.opts)")
#paramtest = IndGenExpansionParams(h,hj,Xtest1,Ptest1,pmidtest1,imp_opt.opts)
#@time IndGenExpansionParams(h,hj,Xtest1,Ptest1,pmidtest1,imp_opt.opts)
#@time IndGenExpansionParams(h,hj,Xtest1,Ptest1,pmidtest1,imp_opt.opts)
println("finished param test")

Profile.init(n = 10^7, delay = 0.003)
@profile IndGenExpansionParams(h,hj,Xtest1,Ptest1,pmidtest1,imp_opt.opts)
Profile.print()

# SETS UP SOLVER OBJECT FOR IMPLICIT FUNCTION
s1 = EAGO_NLPSolver(ImplicitOpts = imp_opt,
                    UBDsolver = IpoptSolver(print_level=0,
                                            tol = 1E5,
                                            nlp_scaling_method = "none",
                                            honor_original_bounds = "yes",
                                            mu_strategy = "adaptive",
                                            mu_oracle = "loqo",
                                            required_infeasibility_reduction = 0.5),
                    probe_depth = -1,
                    variable_depth = 1000,#1000,
                    STD_RR_depth = 10,
                    DAG_depth = -1,
                    atol = 1E-1,
                    rtol = 1E-2)

mKinetic = Model(solver = s1)

z = @variable(mKinetic, [i=1:(5*n+3)], lowerbound=LBDKinetic(i,n), upperbound=UBDKinetic(i,n))

@NLparameter(mKinetic, delT == 0.01)
@NLparameter(mKinetic, K2 == 15339.173673202893)
@NLparameter(mKinetic, K3 == 2*15339.173673202893)
@NLparameter(mKinetic, k1 == 53.0)
@NLparameter(mKinetic, k5 == 1.2E-3)
@NLparameter(mKinetic, k1s == 53.0*1E-6)
@NLparameter(mKinetic, cO2 == 2.0E-3)

@NLconstraint(mKinetic, -z[1] + delT*(k1*z[4]*z[5] - cO2*(z[(5*n+1)]+z[5*n+2])*z[1] + (z[5*n+1]/K2)*z[3] + (z[5*n+2]/K3)*z[2] - k5*z[1]^2) == 0)
@NLconstraint(mKinetic, -z[2] + delT*(z[5*n+2]*cO2*z[1]-(z[5*n+2]/K3+z[(5*n+3)])*z[2]) == 0)
@NLconstraint(mKinetic, -z[3] + delT*(z[5*n+1]*cO2*z[1]-(z[5*n+1]/K2)*z[3]) == 0)
@NLconstraint(mKinetic, 0.4 - z[4] + delT*(-k1*z[4]*z[5]) == 0)
@NLconstraint(mKinetic, 140.0 - z[5] + delT*(-k1*z[4]*z[5]) == 0)

@NLconstraint(mKinetic, [i=2:n], z[5*(i-1)-4] - z[5*(i-1)+1] + delT*(k1*z[5*(i-1)+4]*z[5*(i-1)+5] - cO2*(z[5*n+1]+z[5*n+2])*z[5*(i-1)+1] + (z[5*n+1]/K2)*z[5*(i-1)+3] + (z[5*n+2]/K3)*z[5*(i-1)+2] - k5*z[5*(i-1)+1]^2) == 0)
@NLconstraint(mKinetic, [i=2:n], z[5*(i-1)-3] - z[5*(i-1)+2] + delT*(z[5*n+2]*cO2*z[5*(i-1)+1]-(z[5*n+2]/K3+z[(5*n+3)])*z[5*(i-1)+2]) == 0)
@NLconstraint(mKinetic, [i=2:n], z[5*(i-1)-2] - z[5*(i-1)+3] + delT*(z[5*n+1]*cO2*z[5*(i-1)+1]-(z[5*n+1]/K2)*z[5*(i-1)+3]) == 0)
@NLconstraint(mKinetic, [i=2:n], z[5*(i-1)-1] - z[5*(i-1)+4] + delT*(-k1*z[5*(i-1)+4]*z[5*(i-1)+5]) == 0)
@NLconstraint(mKinetic, [i=2:n], z[5*(i-1)] - z[5*(i-1)+5] + delT*(-k1*z[5*(i-1)+4]*z[5*(i-1)+5]) == 0)

I_Data = [66.0952, 104.762, 110.333, 114.905, 122.238, 125.429, 125.429, 123.476, 121.286,
          118.857, 117.667, 116.143, 113.857, 111.571, 108.81,  105.952, 104.048, 102.048,
          100.143, 98.5238, 96.2381, 94.381,  91.6667, 89.5714, 87.1429, 84.8571, 83.4286,
          81.1905, 78.9048, 77.0476, 75.4762, 73.4762, 71.8095, 70.6667, 68.381,  67.3333,
          65.0952, 63.7143, 62.0476, 60.8571, 59.619,  58.2857, 57.4762, 56.4762, 55.8095,
          54.5238, 53,      51.8571, 50.4286, 49.381,  47.9524, 47.3714, 46.8952, 46.4857,
          45.9048, 45.0762, 44.3238, 43.4143, 43.5429, 42.3619, 41.8381, 40.2381, 39.1286,
          38.7857, 37.081,  36.9524, 36.581,  36.281,  35.3476, 34.8905, 34.1667, 33.6714,
          32.9667, 31.8429, 31.5429, 31.1476, 30.9905, 29.9571, 29.1333, 28.7857, 28.4429,
          28.3476, 27.5429, 27.4333, 27.6048, 27.1762, 27.2,    26.4333, 25.7619, 24.8095,
          24.7429, 24.2857, 24.1714, 23.5667, 23.5476, 23.3952, 22.919,  22.3095, 21.8048,
          21.2857, 21.2048, 20.8429, 20.4429, 20.0048, 19.9381, 19.5,    19.8667, 18.9333,
          19.1381, 18.9619, 18.5476, 17.9048, 17.7571, 18.5333, 18.3762, 18.3571, 18.3286,
          18.2762, 18.3952, 17.5952, 18.1524, 18.1952, 17.8476, 17.9095, 17.5048, 17.5,
          15.9619, 16.2095, 16.181,  15.6952, 15.7095, 15.4619, 15.9476, 16,      16.1952,
          16.1143, 15.7429, 15.5762, 15.7048, 15.8095, 15.6667, 14.9048, 14.5857, 14.7524,
          14.7571, 14.9762, 14.5333, 14.5524, 14.0143, 13.6286, 13.4429, 13.4667, 13.319,
          12.9333, 13.1238, 12.7476, 12.9333, 13.0714, 13.0714, 12.7619, 12.4238, 12.5143,
          12.9143, 12.5714, 13.3667, 13.2286, 13.7905, 13.7571, 13.5905, 12.9667, 12.981,
          12.8857, 12.919,  13.0143, 13.0095, 12.3857, 12.5571, 12.3429, 12.7571, 12.681,
          12.5429, 12.1857, 12.7905, 12.5571, 12.8429, 12.5476, 12.5714, 12.3762, 11.9952,
          11.4571, 11.3,    11.1524, 11.681,  11.619,  11.9048, 12,      12.0762, 11.9143,
          11.7619, 11.5333]

@NLobjective(mKinetic, :Min, sum((z[5*(i-1)+1] + (2.0/21.0)*z[5*(i-1)+2] + (2.0*21.0)*z[5*(i-1)+3] - I_Data[i])^2 for i=1:n))

solve(mKinetic)
