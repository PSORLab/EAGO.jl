module ImplicitBnd

using Compat
using Compat.Test
using StaticArrays
using IntervalArithmetic
using EAGO


function XaboutY(x,y,tol)
    return abs(x-y) <= tol
end

@testset "Test Base Implicit Routines" begin

opts1 =  mc_opts{Float64}(0.5,1,:Dense,:Newton,1,1,1E-10)

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

np = 1
szero = @SVector zeros(np)
sone = @SVector ones(np)
p_mc = [SMCg{np,Interval{Float64},Float64}(p[i],p[i],sone,sone,@interval(P[i].lo,P[i].hi),false,xIBox,mBox) for i=1:np]
param = GenExpansionParams(h1,hj1,X,P,pmid,opts1)

hbnds = MC_impRelax(h1,hj1,p_mc,pmid,X,P,opts1,param)
fbnds = impRelax_f(f,h1,hj1,X,P,p,pmid,opts1,param)
fgbnds = impRelax_fg(f,g,h1,hj1,X,P,p,pmid,opts1,param)
#=
param = GenExpansionParams(h1,hj1,X,P,pmid,opts1)
println("param[2][1].cc: $(param[2][1].cc)")
println("param[2][1].cv: $(param[2][1].cv)")
@test XaboutY(param[2][1].cc,-0.46623,1E-4)
@test XaboutY(param[2][1].cv,-0.57927,1E-4)

hbnds = MC_impRelax(h1,hj1,p_mc,pmid,X,P,opts1,param)
@test XaboutY(hbnds[1].cc,-0.49073,1E-4)
@test XaboutY(hbnds[1].cv,-0.57581,1E-4)

fbnds = impRelax_f(f,h1,hj1,X,P,p,pmid,opts1,param)
@test XaboutY(fbnds.cc,3.95557,1E-4)
@test XaboutY(fbnds.cv,2.99472,1E-4)

fgbnds = impRelax_fg(f,g,h1,hj1,X,P,p,pmid,opts1,param)
@test XaboutY(fgbnds[1].cc,3.95557,1E-4)
@test XaboutY(fgbnds[1].cv,2.99472,1E-4)
@test XaboutY(fgbnds[2][1].cc,3.95557,1E-4)
@test XaboutY(fgbnds[2][1].cv,2.99472,1E-4)
@test XaboutY(fgbnds[2][2].cc ,11.45557,1E-4)
@test XaboutY(fgbnds[2][2].cv,10.49472,1E-4)
=#
end
end
