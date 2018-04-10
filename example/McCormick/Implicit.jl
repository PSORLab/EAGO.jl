#workspace()

using EAGO
using IntervalArithmetic
using StaticArrays

# sets up initial options
opts1 = mc_opts(Float64)          # sets options for relaxtion
opts1.kmax = 3             # sets number of iterations
#opts.style = "NewtonGS"   # sets style of contractor
opts1.style = "KrawczykCW"   # sets style of contractor

#=
generates the expansion point parameters for the function using the opts
options using inverse preconditioner
=#
f(x,p) = x[1]*p[1]+p[1]
g(x,p) = [x[1]*p[1]+p[1];
          x[1]*p[1]+2*p[1]]
function h1(x,p)
    t1 = x[1]^2
    t2 = x[1]*p[1]
    t3 = 4.0
    t4 = t1 + t2
    t5 = t2 + t3
    return [t5]
end
hj1(x,p) = [2*x[1]+p[1]]
P = [Interval(6.0,9.0)]
X = [Interval(-0.78,-0.4)]
p = [7.5]
pmid = mid.(P)

xIntv1 = Interval(1.0,3.0)
xIBox = SVector{1,Interval{Float64}}([xIntv1])
mBox = mid.(xIBox)

#=
relaxes the equality h(x,p)
=#
np = 1
szero = @SVector zeros(np)
sone = @SVector ones(np)
p_mc = [SMCg{np,Interval{Float64},Float64}(p[i],p[i],sone,sone,@interval(P[i].lo,P[i].hi),false,xIBox,mBox) for i=1:np]

#=
generates the expansion point parameters for the function using the opts
options using in place LDU full pivot conditioner with sparse calc
=#
f1(x,p) = x[1]*p[1]+p[1]
g1(x,p) = [x[1]*p[1]+p[1];
          x[1]*p[1]+2*p[1]]
function h1!(hout,x,p)
    t1 = x[1]^2
    t2 = x[1]*p[1]
    t3 = 4.0
    t4 = t1 + t2
    t5 = t2 + t3
    hout[1] = t5
end
function hj1!(hout,x,p)
    hout[1] = 2*x[1]+p[1]
end
P1 = [Interval(6.0,9.0)]
X1 = [Interval(-0.78,-0.4)]
p1 = [7.5]
pmid1 = mid.(P1)

np1 = 1
szero1 = @SVector zeros(np1)
sone1 = @SVector ones(np1)
p_mc1 = [SMCg{np,Interval{Float64},Float64}(p1[i],p1[i],sone1,sone1,@interval(P1[i].lo,P1[i].hi),false,xIBox,mBox) for i=1:np1]

param = GenExpansionParams(h1,hj1,X,P,pmid,opts1)
# param[2][1].cc = -0.46623417721518967
# param[2][1].cv = -0.5792721518987343
hbnds = MC_impRelax(h1,hj1,p_mc,pmid,X,P,opts1,param)
# hbnds[1].cc = -0.4907369110006179
# hbnds[1].cv = -0.5758134311908059
fbnds = impRelax_f(f,h1,hj1,X,P,p,pmid,opts1,param)
# fbnds.cc = 3.955578533996293
# fbnds.cv = 2.9947278349146904
fgbnds = impRelax_fg(f,g,h1,hj1,X,P,p,pmid,opts1,param)
# fgbnds[1].cc = 3.955578533996293
# fgbnds[1].cv = 2.9947278349146904
# fgbnds[2][1].cc = 3.955578533996293
# fgbnds[2][1].cv = 2.9947278349146904
# fgbnds[2][2].cc = 11.455578533996293
# fgbnds[2][2].cv = 10.494727834914691

param1 = InGenExpansionParams(h1!,hj1!,X1,P1,pmid1,opts1)
# param1[2][1].cc = -0.46623417721518967
# param1[2][1].cv = -0.5792721518987343
# relaxation of f(x,p) & g(x,p) at (x(p),p) SPARSE SOLVER
hbnds1 = MC_NimpRelax(h1!,hj1!,p_mc1,pmid1,X1,P1,opts1,param1)
# hbnds1[1].cc = -0.4907369110006179
# hbnds1[1].cv = -0.5758134311908059
fbnds1 = NimpRelax_f(f1,h1!,hj1!,X1,P1,p,pmid1,opts1,param1)
# fbnds1.cc = 3.955578533996293
# fbnds1.cv = 2.9947278349146904
fgbnds1 = NimpRelax_fg(f1,g1,h1!,hj1!,X1,P1,p,pmid1,opts1,param1)
# fgbnds1[1].cc = 3.955578533996293
# fgbnds1[1].cv = 2.9947278349146904
# fgbnds1[2][1].cc = 3.955578533996293
# fgbnds1[2][1].cv = 2.9947278349146904
# fgbnds1[2][2].cc = 11.455578533996293
# fgbnds1[2][2].cv = 10.494727834914691

# relaxation of f(x,p) & g(x,p) at (x(p),p) DIRECT SOLVER
param2 = IndGenExpansionParams(h1!,hj1!,X1,P1,pmid1,opts1)
# param2[2][1].cc = -0.46623417721518967
# param2[2][1].cv = -0.5792721518987343
hbnds2 = MC_NdimpRelax(h1!,hj1!,p_mc1,pmid1,X1,P1,opts1,param2)
# hbnds2[1].cc = -0.4907369110006179
# hbnds2[1].cv = -0.5758134311908059
fbnds2 = NdimpRelax_f(f1,h1!,hj1!,X1,P1,p,pmid1,opts1,param2)
# fbnds2.cc = 3.955578533996293
# fbnds2.cv = 2.9947278349146904
fgbnds2 = NdimpRelax_fg(f1,g1,h1!,hj1!,X1,P1,p,pmid1,opts1,param2)
# fgbnds2[1].cc = 3.955578533996293
# fgbnds2[1].cv = 2.9947278349146904
# fgbnds2[2][1].cc = 3.955578533996293
# fgbnds2[2][1].cv = 2.9947278349146904
# fgbnds2[2][2].cc = 11.455578533996293
# fgbnds2[2][2].cv = 10.494727834914691

opts2 = mc_opts(Float64)   # sets options for relaxtion
opts2.kmax = 3             # sets number of iterations
opts2.style = "NewtonGS"   # sets style of contractor

param3 = GenExpansionParams(h1,hj1,X,P,pmid,opts2)
# param3[2][1].cc = -0.4593323920935859
# param3[2][1].cv = -0.5667248956273347
hbnds3 = MC_impRelax(h1,hj1,p_mc,pmid,X,P,opts2,param)
# hbnds3[1].cc = -0.5014435013618904
# hbnds3[1].cv = -0.5667248956273347
fbnds3 = impRelax_f(f,h1,hj1,X,P,p,pmid,opts2,param)
# fbnds3.cc = 3.891338991828657
# fbnds3.cv = 3.0592452208305847
fgbnds3 = impRelax_fg(f,g,h1,hj1,X,P,p,pmid,opts2,param)
# fgbnds3[1].cc = 3.891338991828657
# fgbnds3[1].cv = 3.0592452208305847
# fgbnds3[2][1].cc = 3.891338991828657
# fgbnds3[2][1].cv = 3.0592452208305847
# fgbnds3[2][2].cc = 11.391338991828658
# fgbnds3[2][2].cv = 10.559245220830585

param4 = InGenExpansionParams(h1!,hj1!,X1,P1,pmid1,opts2)
# param4[2][1].cc = -0.4593323920935859
# param4[2][1].cv = -0.5667248956273347
hbnds4 = MC_NimpRelax(h1!,hj1!,p_mc1,pmid1,X1,P1,opts2,param1)
# hbnds4[1].cc = -0.5014435013618905
# hbnds4[1].cv = -0.5667248956273346
fbnds4 = NimpRelax_f(f1,h1!,hj1!,X1,P1,p,pmid1,opts2,param1)
# fbnds4.cc = 3.891338991828657
# fbnds4.cv = 3.0592452208305865
fgbnds4 = NimpRelax_fg(f1,g1,h1!,hj1!,X1,P1,p,pmid1,opts2,param1)
# fgbnds4[1].cc = 3.891338991828657
# fgbnds4[1].cv = 3.0592452208305865
# fgbnds4[2][1].cc = 3.891338991828657
# fgbnds4[2][1].cv = 3.0592452208305865
# fgbnds4[2][2].cc = 11.391338991828658
# fgbnds4[2][2].cv = 10.559245220830586

# relaxation of f(x,p) & g(x,p) at (x(p),p) DIRECT SOLVER
param5 = IndGenExpansionParams(h1!,hj1!,X1,P1,pmid1,opts2)
# param5[2][1].cc = -0.4593323920935859
# param5[2][1].cv = -0.5667248956273347
hbnds5 = MC_NdimpRelax(h1!,hj1!,p_mc1,pmid1,X1,P1,opts2,param2)
# hbnds5[1].cc = -0.5014435013618904
# hbnds5[1].cv = -0.5667248956273347
fbnds5 = NdimpRelax_f(f1,h1!,hj1!,X1,P1,p,pmid1,opts2,param2)
# fbnds5.cc = 3.891338991828657
# fbnds5.cv = 3.0592452208305847
fgbnds5 = NdimpRelax_fg(f1,g1,h1!,hj1!,X1,P1,p,pmid1,opts2,param2)
# fgbnds5[1].cc = 3.891338991828657
# fgbnds5[1].cv = 3.0592452208305847
# fgbnds5[2][1].cc = 3.891338991828657
# fgbnds5[2][1].cv = 3.0592452208305847
# fgbnds5[2][2].cc = 11.391338991828658
# fgbnds5[2][2].cv = 10.559245220830585
