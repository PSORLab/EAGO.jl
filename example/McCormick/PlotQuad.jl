workspace()
using EAGO
using IntervalArithmetic
using StaticArrays
#using DataFrames
#using CSV
using Plots
#using PyPlot
pyplot()

xl = 6.0
xu = 9.0

intv = 0.05
nx = Integer((xu-xl)/intv+1)
opts1 =  mc_opts{Float64}(0.5,1,:Dense,:Newton,1,1,0.0)
function h1(x,p)
    t1 = x[1]^2
    t2 = x[1]*p[1]
    t3 = 4.0
    t4 = t1 + t2
    t5 = t2 + t3
    return [t5]
end
hj1(x,p) = [2.0*x[1]+p[1]]
P = [Interval(6.0,9.0)]
pmid = [7.5]
X = [Interval(-0.78,-0.4)]

x_grid = zeros(Float64,nx)
values = zeros(Float64,nx)
solbranch = zeros(Float64,nx)

ccvalNSnr = zeros(Float64,nx)  # Storage for nonsmooth, subgradient refinement
cvvalNSnr = zeros(Float64,nx)
ccvalNSnn = zeros(Float64,nx)  # Storage for nonsmooth, no subgradient refinement
cvvalNSnn = zeros(Float64,nx)

ccvalNSsr = zeros(Float64,nx)  # Storage for smooth, subgradient refinement
cvvalNSsr = zeros(Float64,nx)
ccvalNSsn = zeros(Float64,nx)  # Storage for smooth, no subgradient refinement
cvvalNSsn = zeros(Float64,nx)

tempbox = SVector{1,Interval{Float64}}(1.0*[Interval(xl,xu)])
tempmid =  mid.(tempbox)

EAGO.set_diff_relax(0)
EAGO.set_subgrad_refine(false)
param_ns_sr = GenExpansionParams(h1,hj1,X,P,pmid,opts1)
temp1 = SMCg{1,Interval{Float64},Float64}[SMCg{1,Interval{Float64},Float64}(Float64(xl),Float64(xl),seed_g(Float64,1,1),
                         seed_g(Float64,1,1),Interval(xl,xu), false, tempbox,tempmid)]
temp1_MC = MC_impRelax(h1,hj1,temp1,pmid,X,P,opts1,param_ns_sr)
count1 = 1
for j=xl:intv:xu
    values[count1] = j
    solbranch[count1] = (-j+sqrt(j^2-16.0))/2.0
    temp1 = SMCg{1,Interval{Float64},Float64}[SMCg{1,Interval{Float64},Float64}(Float64(j),Float64(j),seed_g(Float64,1,1),
                             seed_g(Float64,1,1),Interval(xl,xu),
                             false, tempbox,tempmid)]
    temp1_MC = MC_impRelax(h1,hj1,temp1,pmid,X,P,opts1,param_ns_sr)
    ccvalNSnn[count1] = temp1_MC[end].cc
    cvvalNSnn[count1] = temp1_MC[end].cv
    count1 += 1
end

EAGO.set_subgrad_refine(false)
count1 = 1
param_ns_nr = GenExpansionParams(h1,hj1,X,P,pmid,opts1)
for j=xl:intv:xu
    temp1 = SMCg{1,Interval{Float64},Float64}[SMCg{1,Interval{Float64},Float64}(Float64(j),Float64(j),seed_g(Float64,1,1),
                                                  seed_g(Float64,1,1),Interval(xl,xu),
                                                  false, tempbox,tempmid)]
    temp1_MC = MC_impRelax(h1,hj1,temp1,pmid,X,P,opts1,param_ns_sr)
    ccvalNSnr[count1] = temp1_MC[end].cc
    cvvalNSnr[count1] = temp1_MC[end].cv
    count1 += 1
end

EAGO.set_diff_relax(1)
EAGO.set_subgrad_refine(true)
param_diff_nr = GenExpansionParams(h1,hj1,X,P,pmid,opts1)
count1 = 1
for j=xl:intv:xu
    temp1 = SMCg{1,Interval{Float64},Float64}[SMCg{1,Interval{Float64},Float64}(Float64(j),Float64(j),seed_g(Float64,1,1),
                             seed_g(Float64,1,1),Interval(xl,xu),
                             false, tempbox,tempmid)]
    temp1_MC = MC_impRelax(h1,hj1,temp1,pmid,X,P,opts1,param_diff_nr)
    ccvalNSsr[count1] = temp1_MC[end].cc
    cvvalNSsr[count1] = temp1_MC[end].cv
    count1 += 1
end
EAGO.set_subgrad_refine(false)
param_diff_sr = GenExpansionParams(h1,hj1,X,P,pmid,opts1)
count1 = 1
for j=xl:intv:xu
    temp1 = SMCg{1,Interval{Float64},Float64}[SMCg{1,Interval{Float64},Float64}(Float64(j),Float64(j),seed_g(Float64,1,1),
                             seed_g(Float64,1,1),Interval(xl,xu),
                             false, tempbox,tempmid)]
    temp1_MC = MC_impRelax(h1,hj1,temp1,pmid,X,P,opts1,param_diff_sr)
    ccvalNSsn[count1] = temp1_MC[end].cc
    cvvalNSsn[count1] = temp1_MC[end].cv
    count1 += 1
end

# TO DO:
# Debug IMPLICIT
# Plotting Code
# Subgradient refinement code (mult,power)

# plots affine relaxations
#=
plot(values,[ccaffNSnr,ccaffNSnn,ccaffNSsr,ccaffNSsn,
             cvaffNSnr,cvaffNSnn,cvaffNSsr,cvaffNSsn,
             solbranch,solbranch,solbranch,solbranch],layout=4)
=#

# plots functions
#plot(Plots.fakedata(50,5),w=3)

plot(values,[ccvalNSnr,ccvalNSnn,ccvalNSsr,ccvalNSsn,
             cvvalNSnr,cvvalNSnn,cvvalNSsr,cvvalNSsn,
             solbranch,solbranch,solbranch,solbranch],layout=4)

gui()
