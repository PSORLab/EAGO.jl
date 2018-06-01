workspace()

using EAGO
using IntervalArithmetic
using StaticArrays

EAGO.set_diff_relax(0)
opts1 = mc_opts{Float64}(0.5,1,:Dense,:Newton,1,1,1E-10)
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
    t5 = t4 + t3
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
nxi = 1
p_mc = [HybridMC{np,Interval{Float64},Float64}(SMCg{np,Interval{Float64},Float64}(p[i],p[i],sone,sone,Interval{Float64}(P[i].lo,P[i].hi),false)) for i=1:np]
xa_mc = HybridMC{np,Interval{Float64},Float64}[HybridMC{np,Interval{Float64},Float64}(SMCg{np,Interval{Float64},Float64}(X[i].lo,X[i].lo,szero,szero,Interval{Float64}(X[i].lo,X[i].lo),false)) for i=1:nxi]
xA_mc = HybridMC{np,Interval{Float64},Float64}[HybridMC{np,Interval{Float64},Float64}(SMCg{np,Interval{Float64},Float64}(X[i].hi,X[i].hi,szero,szero,Interval{Float64}(X[i].hi,X[i].hi),false)) for i=1:nxi]
z_mc = 0.5*xa_mc+(one(Float64)-0.5)*xA_mc
h_mc = h1(z_mc,p_mc)

param = GenExpansionParams(h1,hj1,X,P,pmid,opts1)

hbnds = MC_impRelax(h1,hj1,p_mc,pmid,X,P,opts1,param)
fbnds = impRelax_f(f,h1,hj1,X,P,p,pmid,opts1,param)
fgbnds = impRelax_fg(f,g,h1,hj1,X,P,p,pmid,opts1,param)

param = GenExpansionParams(h1,hj1,X,P,pmid,opts1)
# param[2][1].cc = -0.46623417721518967
# param[2][1].cv = -0.5792721518987343
#println("typeof param_ns_sr: $(typeof(param_ns_sr))")
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

opts2 =  mc_opts{Float64}(0.5,3,:Dense,:Newton,1,1,0.0)

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

#=
np = 2
opts3 =  mc_opts{Float64}(0.5,3,:Dense,:Krawczyk,2,2,0.0)
opts4 =  mc_opts{Float64}(0.5,3,:Dense,:Newton,2,2,0.0)
P2 = [Interval(5.0,7.0),Interval(5.0,7.0)]
pmid2 = mid.(P2)
p2 = copy(pmid2)
Z2 = [Interval(-1.5, 0.0),Interval(0.0, 0.5)]
Z2c = copy(Z2)
h2(z,p) = [z[1]^2+z[2]^2+p[1]*z[1]+4.0;
           z[1]+p[2]*z[2]]
hj2(z,p) = [(z[1]+p[1]) (z[2]);
              one(p[1])  p[2]]
sone = @SVector ones(np)
xIntv1 = Interval(1.0,3.0)
xIBox2 = SVector{2,Interval{Float64}}([xIntv1, xIntv1])
mBox2 = mid.(xIBox2)
p_mc2 = [SMCg{np,Interval{Float64},Float64}(p2[i],p2[i],sone,sone,@interval(P2[i].lo,P2[i].hi),false,xIBox2,mBox2) for i=1:np]
param4 = GenExpansionParams(h2,hj2,Z2,P2,pmid2,opts3)
hbnds4 = MC_impRelax(h2,hj2,p_mc2,pmid2,Z2,P2,opts3,param4)
fbnds4 = impRelax_f(f,h2,hj2,Z2,P2,p2,pmid2,opts3,param4)
fgbnds4 = impRelax_fg(f,g,h2,hj2,Z2,P2,p2,pmid2,opts3,param4)
param5 = GenExpansionParams(h2,hj2,Z2c,P2,pmid2,opts4)
hbnds5 = MC_impRelax(h2,hj2,p_mc2,pmid2,Z2c,P2,opts4,param5)
fbnds5 = impRelax_f(f,h2,hj2,Z2c,P2,p2,pmid2,opts4,param5)
fgbnds5 = impRelax_fg(f,g,h2,hj2,Z2c,P2,p2,pmid2,opts4,param5)
=#
