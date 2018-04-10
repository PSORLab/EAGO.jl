#workspace()
#=
using EAGOParametricInterval
using IntervalArithmetic

opt1 = Param_Bisect_Opts() # sets bisection options to default
opt1.kmax_main = 1      # number of iterations for generalized bisection routine
opt1.kmax_cntr = 100       # number of iterations for contractor
opt1.style = "NewtonGS"    # specifies NewtonGS
opt1.DAGflag = false       # disable interval constraint prop

opt2 = Param_Bisect_Opts() # sets bisection options to default
opt2.kmax_main = 1      # number of iterations for generalized bisection routine
opt2.kmax_cntr = 100        # number of iterations for contractor
opt2.style = "KrawczykCW"    # specifies NewtonGS
opt2.DAGflag = false       # disable interval constraint prop

h1(z,p) = [z[1]^2+z[2]^2+p[1]*z[1]+4.0;
           z[1]+p[2]*z[2]]
h1j(z,p) = [2.0*z[1]+p[1]   2.0*z[2];
                      1.0   p[2]]
X1 = [Interval(-10.0,10.0),Interval(-2.0,2.0)]
P1 = [Interval(5.0,7.0),Interval(5.0,7.0)]

println("-------------------------------------------------")
println("----------- EXAMPLE 1 SUMMARY (NEWTON) ----------")
println("-------------------------------------------------")
println(" ")
Nsol1,Nin1,Niter1,Nhist1 = Generalized_Param_Bisection(X1,P1,h1,h1j,opt1)

println("-------------------------------------------------")
println("----------- EXAMPLE 1 SUMMARY (KRAWCZ) ----------")
println("-------------------------------------------------")
println(" ")
Ksol1,Kin1,Kiter1,Khist1 = Generalized_Param_Bisection(X1,P1,h1,h1j,opt2)

println("----------- EXAMPLE 2 ----------")
h2(z,p) = [(3.25-z[1])/p[1]-z[3];
           z[1]/p[2]-z[3];
           z[2]-(z[1]^2)/(1+z[1]^2)]
h2j(z,p) = [-1/p[1] 0 (-1);
            1/p[2] 0 (-1);
            -2*z[1]/(1+z[1]^2)^2 1 0]
X2 = [Interval(-30..30),Interval(-30..30),Interval(-30..30)]
P2 = [Interval(1800..2200),Interval(900..1100)]
# #=
println("-------------------------------------------------")
println("----------- EXAMPLE 2 SUMMARY (NEWTON) ----------")
println("-------------------------------------------------")
println(" ")
#Nsol2,Nin2,Niter2,Nhist2 = Generalized_Param_Bisection(X2,P2,h2,h2j,opt1)

println("-------------------------------------------------")
println("----------- EXAMPLE 2 SUMMARY (KRAWCZ) ----------")
println("-------------------------------------------------")
println(" ")
#Ksol2,Kin2,Kiter2,Khist2 = Generalized_Param_Bisection(X2,P2,h2,h2j,opt2)


h3(z,p) = [(exp(38.0*z[1]-1.0+log(10.0^(-9.0)))) + p[1]*z[1]-1.6722*z[2] + 0.6689*z[3]-8.0267
           (exp(38.0*z[2]-1.0)+log(1.98*10^(-9.0))) + 0.6622*z[1] + p[2]*z[2] + 0.6622*z[3] + 4.0535
           (exp(38.0*z[3]- 1.0+(10.0^(-9.0)))) + z[1] - z[2] + p[3]*z[3] - 6.0]
h3j(z,p) = [p[1]+exp(38.0*z[1]-1.0+log(38.0*(10.0^(-9.0)))) -1.6722  0.6689;
            0.6622 exp(38.0*z[2]-1.0+log(38*(1.98*10.0^(-9.0))))+p[2] 0.6622;
            1.0 -2.0 p[3]+exp(38.0*z[3]-1.0+log(38.0*(10.0^(-9.0))));]
X3 = [Interval(-30.0,30.0),Interval(-30.0,30.0),Interval(-30.0,30.0)]
P3 = [Interval(0.6020,0.7358),Interval(1.2110,1.4801),Interval(3.6,4.4)]
h3temp = h3(X3,P3)
hj3temp = h3j(X3,P3)
Nsol3,Nin3,Niter3,Nhist3 = Generalized_Param_Bisection(X3,P3,h3,h3j,opt1)
Ksol3,Kin3,Kiter3,Khist3 = Generalized_Param_Bisection(X3,P3,h3,h3j,opt2)

h4(z,p) = [(-z[1]^2+p[1])*z[1]]
h4j(z,p) = [-3*z[1]^2+p[1]]
X4 = [Interval(-10.0,10.0)]
P4 = [Interval(0.25..20)]
#Nsol4,Nin4,Niter4,Nhist4 = Generalized_Param_Bisection(X4,P4,h4,h4j,opt1)
#Ksol4,Kin4,Kiter4,Khist4 = Generalized_Param_Bisection(X4,P4,h4,h4j,opt2)

h5(z,p) = [p[1]*z[1]^5-25.2*z[1]^3+6*p[1]*z[1]-p[2]*z[2];
           2*p[2]*z[2]-p[1]*z[1]]
h5j(z,p) = [5*p[1]*z[1]^4-3*25.2*z[1]^2+6*p[1]  (-p[2]);
            -p[1]  2*p[2]]
X5 = [Interval(3.0,5.0),Interval(4.25,7.75)]
P5 = [Interval(-10.0,10.0),Interval(-10.0,10.0)]
#Nsol5,Nin5,Niter5,Nhist5 = Generalized_Param_Bisection(X5,P5,h5,h5j,opt1)
#Ksol5,Kin5,Kiter5,Khist5 = Generalized_Param_Bisection(X5,P5,h5,h5j,opt2)



println("---------- BEGIN SUMMARY OF TEST PROBLEMS --------------")


println("Test Problem #1: ")
println("Newton   Iterations Taken: $(Niter1) (Target = 169)")
println("Krawczyk Iterations Taken: $(Kiter1) (Target = 131)")
println("Newton   Boxes Taken: $(length(Nsol1)) (Target = 5)")
println("Krawczyk Boxes Taken: $(length(Ksol1)) (Target = 3)")

#=
println("Test Problem #2: ")
println("Newton   Iterations Taken: $(Niter2) (Target = 1)")
println("Krawczyk Iterations Taken: $(Kiter2) (Target = 1)")
println("Newton   Boxes Taken: $(length(Nsol2)) (Target = 1)")
println("Krawczyk Boxes Taken: $(length(Ksol2)) (Target = 1)")
=#
#=
println("Test Problem #3: ")
println("Newton   Iterations Taken: $(Niter3) (Target = 55)")
println("Krawczyk Iterations Taken: $(Kiter3) (Target = 329)")
println("Newton   Boxes Taken: $(length(Nsol3)) (Target = 1)")
println("Krawczyk Boxes Taken: $(length(Ksol3)) (Target = 1)")

println("Test Problem #4: ")
println("Newton   Iterations Taken: $(Niter4) (Target = 270)")
println("Krawczyk Iterations Taken: $(Kiter4) (Target = 295)")
println("Newton   Boxes Taken: $(length(Nsol4)) (Target = 15)")
println("Krawczyk Boxes Taken: $(length(Ksol4)) (Target = 15)")


println("Test Problem #5: ")
println("Newton   Iterations Taken: $(Niter5) (Target = 2392)")
println("Krawczyk Iterations Taken: $(Kiter5) (Target = 2599)")
println("Newton   Boxes Taken: $(length(Nsol5)) (Target = 80)")
println("Krawczyk Boxes Taken: $(length(Ksol5)) (Target = 80)")


println("---------- END SUMMARY OF TEST PROBLEMS --------------")


# generates directed graph contractor parameters
if (DAGflag)
  DAGr::Int64 = opt.DAGpass
  DAGpack::Vector{Interval{Float64}} = vcat(Xin,Pin)
  exprs::Vector{Expr} = vcat(opt.DAGh,opt.DAGg)
  hbnds::Vector{Float64} = zeros(Float64,length(Xin))
  if (length(opt.DAGgL)>0)
   gL::Vector{Float64} = vcat(hbnds,opt.DAGgL)
   gU::Vector{Float64} = vcat(hbnds,opt.DAGgU)
  else
   gL = vcat(hbnds)
   gU = vcat(hbnds)
  end
  npx::Int64 = length(Xin) + length(Pin)
  DAGparam = Generate_TapeList(exprs,npx,gL,gU)
end
=#

#=
if (~Eflag)
 if (DAGflag)
    #### Contracts box using interval contractors on DAG ####
    DAGpack = vcat(Xw,Pw)
    DAGContractor!(DAGpack,DAGparam,DAGr)
    for i=1:length(DAGpack)
      #### Fathoms by interval contractors ####
      if (isempty(DAGpack[i]))
       Eflag = true
       break
      end
    end
    for i=1:length(Xw)
      Xw[i] = DAGpack[i]
    end
    count = 1
    for i=(length(Xw)+1):(length(Xw)+length(Pw))
      Pw[count] = DAGpack[i]
      count += 1
    end
 end
end
=#
