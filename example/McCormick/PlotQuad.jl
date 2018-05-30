workspace()
using EAGO
using DataFrames
using IntervalArithmetic
using StaticArrays
#using DataFrames
using CSV
using Plots
#using PyPlot
pyplot()

xl = 6.0
xu = 9.0

yl = -1.00
yu = 1.0

intv = 0.05
nx = Integer((xu-xl)/intv+1)
opts1 =  mc_opts{Float64}(0.5,2,:Dense,:Newton,1,1,0.0)
function h1(x,p)
    return [x[1]^2-x[1]*p[1]+4.0]
end
hj1(x,p) = [2.0*x[1]-p[1]]
function h2(x,p)
    return [x[1]^2+x[1]*p[1]+4.0]
end
hj2(x,p) = [2.0*x[1]+p[1]]

P = [Interval(6.0,9.0)]
pmid = [7.5]
X = [Interval(yl,yu)]

x_grid = zeros(Float64,nx)
values = zeros(Float64,nx)
solbranch = zeros(Float64,nx)
solbranch2 = zeros(Float64,nx)

sto_cc_grad = [0.0]
sto_cv_grad = [0.0]
sto_cc_grad1 = [0.0]
sto_cv_grad1 = [0.0]

sto_cc_grada = [0.0]
sto_cv_grada = [0.0]
sto_cc_grad1a = [0.0]
sto_cv_grad1a = [0.0]

ccvalNSnr = zeros(Float64,nx)  # Storage for nonsmooth, subgradient refinement
cvvalNSnr = zeros(Float64,nx)
ccvalNSnn = zeros(Float64,nx)  # Storage for nonsmooth, no subgradient refinement
cvvalNSnn = zeros(Float64,nx)

ccvalNSsr = zeros(Float64,nx)  # Storage for smooth, subgradient refinement
cvvalNSsr = zeros(Float64,nx)
ccvalNSsn = zeros(Float64,nx)  # Storage for smooth, no subgradient refinement
cvvalNSsn = zeros(Float64,nx)

ccvalNSnr2 = zeros(Float64,nx)  # Storage for nonsmooth, subgradient refinement
cvvalNSnr2= zeros(Float64,nx)
ccvalNSnn2 = zeros(Float64,nx)  # Storage for nonsmooth, no subgradient refinement
cvvalNSnn2 = zeros(Float64,nx)

ccvalNSsr2 = zeros(Float64,nx)  # Storage for smooth, subgradient refinement
cvvalNSsr2 = zeros(Float64,nx)
ccvalNSsn2 = zeros(Float64,nx)  # Storage for smooth, no subgradient refinement
cvvalNSsn2 = zeros(Float64,nx)

tempbox = SVector{1,Interval{Float64}}(1.0*[Interval(xl,xu)])
tempmid =  mid.(tempbox)

EAGO.set_diff_relax(0)
EAGO.set_subgrad_refine(false)
param_ns_sr = GenExpansionParams(h1,hj1,X,P,pmid,opts1)
temp1 = SMCg{1,Interval{Float64},Float64}[SMCg{1,Interval{Float64},Float64}(Float64(xl),Float64(xl),seed_g(Float64,1,1),
                         seed_g(Float64,1,1),Interval(xl,xu), false)]
temp1_MC = MC_impRelax(h1,hj1,temp1,pmid,X,P,opts1,param_ns_sr)
count1 = 1
for j=xl:intv:xu
    values[count1] = j
    solbranch[count1] = (j-sqrt(j^2-16.0))/2.0
    solbranch2[count1] = (-j+sqrt(j^2-16.0))/2.0
    temp1 = SMCg{1,Interval{Float64},Float64}[SMCg{1,Interval{Float64},Float64}(Float64(j),Float64(j),seed_g(Float64,1,1),
                             seed_g(Float64,1,1),Interval(xl,xu),
                             false)]
    temp1_MC = MC_impRelax(h1,hj1,temp1,pmid,X,P,opts1,param_ns_sr)
    ccvalNSnn[count1] = temp1_MC[end].cc
    cvvalNSnn[count1] = temp1_MC[end].cv
    count1 += 1
end

EAGO.set_subgrad_refine(false)
count1 = 1
param_ns_nr1 = GenExpansionParams(h1,hj1,X,P,pmid,opts1)
param_ns_nr2 = GenExpansionParams(h2,hj2,X,P,pmid,opts1)
for j=xl:intv:xu
    temp1 = SMCg{1,Interval{Float64},Float64}[SMCg{1,Interval{Float64},Float64}(Float64(j),Float64(j),seed_g(Float64,1,1),
                                                  seed_g(Float64,1,1),Interval(xl,xu),
                                                  false)]
    temp1_MC = MC_impRelax(h1,hj1,temp1,pmid,X,P,opts1,param_ns_nr1)
    temp2_MC = MC_impRelax(h2,hj2,temp1,pmid,X,P,opts1,param_ns_nr2)
    if (count1 == 31)
        println("temp1 at 31: $(temp1_MC)")
        sto_cc_grad[1] = temp1_MC[1].cc_grad[1]
        sto_cv_grad[1] = temp1_MC[1].cv_grad[1]
        sto_cc_grada[1] = temp2_MC[1].cc_grad[1]
        sto_cv_grada[1] = temp2_MC[1].cv_grad[1]
    end
    ccvalNSnr[count1] = temp1_MC[end].cc
    cvvalNSnr[count1] = temp1_MC[end].cv
    ccvalNSnr2[count1] = temp2_MC[end].cc
    cvvalNSnr2[count1] = temp2_MC[end].cv
    count1 += 1
end

EAGO.set_diff_relax(1)
EAGO.set_subgrad_refine(true)
param_diff_nr1 = GenExpansionParams(h1,hj1,X,P,pmid,opts1)
param_diff_nr2 = GenExpansionParams(h2,hj2,X,P,pmid,opts1)
count1 = 1
for j=xl:intv:xu
    temp1 = SMCg{1,Interval{Float64},Float64}[SMCg{1,Interval{Float64},Float64}(Float64(j),Float64(j),seed_g(Float64,1,1),
                             seed_g(Float64,1,1),Interval(xl,xu),
                             false)]
    temp1_MC = MC_impRelax(h1,hj1,temp1,pmid,X,P,opts1,param_diff_nr1)
    temp2_MC = MC_impRelax(h2,hj2,temp1,pmid,X,P,opts1,param_diff_nr2)
    if (count1 == 31)
        println("temp1 at 31: $(temp1_MC)")
        sto_cc_grad1[1] = temp1_MC[1].cc_grad[1]
        sto_cv_grad1[1] = temp1_MC[1].cv_grad[1]
        sto_cc_grad1a[1] = temp2_MC[1].cc_grad[1]
        sto_cv_grad1a[1] = temp2_MC[1].cv_grad[1]
    end
    ccvalNSsr[count1] = temp1_MC[end].cc
    cvvalNSsr[count1] = temp1_MC[end].cv
    ccvalNSsr2[count1] = temp2_MC[end].cc
    cvvalNSsr2[count1] = temp2_MC[end].cv
    count1 += 1
end
EAGO.set_subgrad_refine(false)
param_diff_sr = GenExpansionParams(h1,hj1,X,P,pmid,opts1)
count1 = 1
for j=xl:intv:xu
    temp1 = SMCg{1,Interval{Float64},Float64}[SMCg{1,Interval{Float64},Float64}(Float64(j),Float64(j),seed_g(Float64,1,1),
                             seed_g(Float64,1,1),Interval(xl,xu),
                             false)]
    temp1_MC = MC_impRelax(h1,hj1,temp1,pmid,X,P,opts1,param_diff_sr)
    ccvalNSsn[count1] = temp1_MC[end].cc
    cvvalNSsn[count1] = temp1_MC[end].cv
    count1 += 1
end


ns_cv_affine = zeros(Float64,nx)
ns_cc_affine = zeros(Float64,nx)
ss_cv_affine = zeros(Float64,nx)
ss_cc_affine = zeros(Float64,nx)
ns_cv_affine_f2 = zeros(Float64,nx)
ns_cc_affine_f2 = zeros(Float64,nx)
ss_cv_affine_f2 = zeros(Float64,nx)
ss_cc_affine_f2 = zeros(Float64,nx)

count1 = 1
for j=xl:intv:xu
    ns_cv_affine[count1] = sto_cv_grad[1]*(j-7.5) + cvvalNSnr[31]
    ns_cc_affine[count1] = sto_cc_grad[1]*(j-7.5) + ccvalNSnr[31]
    ss_cv_affine[count1] = sto_cv_grad1[1]*(j-7.5) + cvvalNSsr[31]
    ss_cc_affine[count1] = sto_cc_grad1[1]*(j-7.5) + ccvalNSsr[31]
    ns_cv_affine_f2[count1] = sto_cv_grada[1]*(j-7.5) + cvvalNSnr2[31]
    ns_cc_affine_f2[count1] = sto_cc_grada[1]*(j-7.5) + ccvalNSnr2[31]
    ss_cv_affine_f2[count1] = sto_cv_grad1a[1]*(j-7.5) + cvvalNSsr2[31]
    ss_cc_affine_f2[count1] = sto_cc_grad1a[1]*(j-7.5) + ccvalNSsr2[31]
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

yl1 = 0.0
yu1 = 1.0
yl2 = -1.0
yu2 = 0.0

plot(values,[#ccvalNSnr,ccvalNSsr,ccvalNSnr2,ccvalNSsr2,
             cvvalNSnr,cvvalNSsr,cvvalNSnr2,cvvalNSsr2,
             ns_cv_affine,ss_cv_affine,ns_cv_affine_f2,ss_cv_affine_f2,
             solbranch,solbranch,solbranch2,solbranch2],layout=4,
             xlims = [(6,9) (6,9) (6,9) (6,9)],
             ylims = [(yl1,yu1) (yl1,yu1) (yl2,yu2) (yl2,yu2)])


#output_frame
#=
[cvvalNSnn,cvvalNSsr,cvvalNSnn2,cvvalNSsr2,
             ns_cv_affine,ss_cv_affine,ns_cv_affine_f2,ss_cv_affine_f2,
             solbranch,solbranch,solbranch2,solbranch2]
=#
#=
plot(values,[solbranch solbranch solbranch2 solbranch2],layout=4,
             xlims = [(6,9) (6,9) 6,9) (6,9)],
             ylims = [(yl1,yu1) (yl1,yu1) (yl2,yu2) (yl2,yu2)])
=#
gui()
