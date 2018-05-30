module Imp_Util_Test

using Compat
using Compat.Test
using IntervalArithmetic
using EAGO
using StaticArrays

# MC Option functions
x = mc_opts()
y = mc_opts(Float64)
x.lambda = 0.75
y.lambda = 0.75
espv = 0.5
seed1a = 0.5*seed_g(Float64,1,2)
intvbox = SVector{2,Interval{Float64}}([Interval(0.5,3.0),Interval(0.5,3.0)])

set_default!(x)
@test x.lambda == 0.5
@test y.lambda == 0.75

# Test Cutting Routines
#Final_Cut(x_mc,x_mc_int)
#Final_Cut(x_mc,x_mc_int)
#set_diff_relax(val::Integer)

# Test Rounding Routines
z_mc = [SMCg{2,Interval{Float64},Float64}(2.0,1.0,seed1a,seed1a,Interval(0.5,3.0),true) for i=1:2]
z_mc1 = EAGO.Rnd_Out_Z_Intv(z_mc,espv)
z_mc2 = [SMCg{2,Interval{Float64},Float64}(2.0,1.0,seed1a,seed1a,Interval{Float64}(0.5-Float64(espv),3.0+Float64(espv)),true) for i=1:2]

@test z_mc1[1].Intv.lo == z_mc2[1].Intv.lo
@test z_mc1[1].Intv.hi == z_mc2[1].Intv.hi
@test z_mc1[2].Intv.lo == z_mc2[2].Intv.lo
@test z_mc1[2].Intv.hi == z_mc2[2].Intv.hi

z_mc = [SMCg{2,Interval{Float64},Float64}(2.0,1.0,seed1a,seed1a,Interval(0.5,3.0),true) for i=1:2]
z_mc1 = EAGO.Rnd_Out_Z_All(z_mc,espv)
z_mc2 = [SMCg{2,Interval{Float64},Float64}(2.0+espv,1.0-espv,seed1a,seed1a,Interval(0.5-espv,3.0+espv),true) for i=1:2]

@test z_mc1[1].Intv.lo == z_mc2[1].Intv.lo
@test z_mc1[1].Intv.hi == z_mc2[1].Intv.hi
@test z_mc1[2].Intv.lo == z_mc2[2].Intv.lo
@test z_mc1[2].Intv.hi == z_mc2[2].Intv.hi
@test z_mc1[1].cc == z_mc2[1].cc
@test z_mc1[1].cv == z_mc2[1].cv
@test z_mc1[2].cc == z_mc2[2].cc
@test z_mc1[2].cv == z_mc2[2].cv

z_mc = [SMCg{2,Interval{Float64},Float64}(2.0,1.0,seed1a,seed1a,Interval(0.5,3.0),true) for i=1:2]
Y_mc = [SMCg{2,Interval{Float64},Float64}(2.0,1.0,seed1a,seed1a,Interval(0.5,3.0),true) for i=1:2,j=1:2]
z_mc1, Y_mc1 = EAGO.Rnd_Out_H_Intv(z_mc,Y_mc,espv)
z_mc2 = [SMCg{2,Interval{Float64},Float64}(2.0,1.0,seed1a,seed1a,Interval(0.5-espv,3.0+espv),true) for i=1:2]
Y_mc2 = [SMCg{2,Interval{Float64},Float64}(2.0,1.0,seed1a,seed1a,Interval(0.5-espv,3.0+espv),true) for i=1:2,j=1:2]

@test z_mc1[1].Intv.lo == z_mc2[1].Intv.lo
@test z_mc1[1].Intv.hi == z_mc2[1].Intv.hi
@test z_mc1[2].Intv.lo == z_mc2[2].Intv.lo
@test z_mc1[2].Intv.hi == z_mc2[2].Intv.hi
@test z_mc1[1].cc == z_mc2[1].cc
@test z_mc1[1].cv == z_mc2[1].cv
@test z_mc1[2].cc == z_mc2[2].cc
@test z_mc1[2].cv == z_mc2[2].cv

@test Y_mc1[1,1].Intv.lo == Y_mc2[1,1].Intv.lo
@test Y_mc1[1,1].Intv.hi == Y_mc2[1,1].Intv.hi
@test Y_mc1[2,1].Intv.lo == Y_mc2[2,1].Intv.lo
@test Y_mc1[2,1].Intv.hi == Y_mc2[2,1].Intv.hi
@test Y_mc1[1,1].cc == Y_mc2[1,1].cc
@test Y_mc1[1,1].cv == Y_mc2[1,1].cv
@test Y_mc1[2,1].cc == Y_mc2[2,1].cc
@test Y_mc1[2,1].cv == Y_mc2[2,1].cv
@test Y_mc1[1,2].Intv.lo == Y_mc2[1,2].Intv.lo
@test Y_mc1[1,2].Intv.hi == Y_mc2[1,2].Intv.hi
@test Y_mc1[2,2].Intv.lo == Y_mc2[2,2].Intv.lo
@test Y_mc1[2,2].Intv.hi == Y_mc2[2,2].Intv.hi
@test Y_mc1[1,2].cc == Y_mc2[1,2].cc
@test Y_mc1[1,2].cv == Y_mc2[1,2].cv
@test Y_mc1[2,2].cc == Y_mc2[2,2].cc
@test Y_mc1[2,2].cv == Y_mc2[2,2].cv

z_mc = [SMCg{2,Interval{Float64},Float64}(2.0,1.0,seed1a,seed1a,Interval(0.5,3.0),true) for i=1:2]
Y_mc = [SMCg{2,Interval{Float64},Float64}(2.0,1.0,seed1a,seed1a,Interval(0.5,3.0),true) for i=1:2,j=1:2]
z_mc1, Y_mc1 = EAGO.Rnd_Out_H_All(z_mc,Y_mc,espv)
z_mc2 = [SMCg{2,Interval{Float64},Float64}(2.0+espv,1.0-espv,seed1a,seed1a,Interval(0.5-espv,3.0+espv),true) for i=1:2]
Y_mc2 = [SMCg{2,Interval{Float64},Float64}(2.0+espv,1.0-espv,seed1a,seed1a,Interval(0.5-espv,3.0+espv),true) for i=1:2,j=1:2]

@test z_mc1[1].Intv.lo == z_mc2[1].Intv.lo
@test z_mc1[1].Intv.hi == z_mc2[1].Intv.hi
@test z_mc1[2].Intv.lo == z_mc2[2].Intv.lo
@test z_mc1[2].Intv.hi == z_mc2[2].Intv.hi
@test z_mc1[1].cc == z_mc2[1].cc
@test z_mc1[1].cv == z_mc2[1].cv
@test z_mc1[2].cc == z_mc2[2].cc
@test z_mc1[2].cv == z_mc2[2].cv

@test Y_mc1[1,1].Intv.lo == Y_mc2[1,1].Intv.lo
@test Y_mc1[1,1].Intv.hi == Y_mc2[1,1].Intv.hi
@test Y_mc1[2,1].Intv.lo == Y_mc2[2,1].Intv.lo
@test Y_mc1[2,1].Intv.hi == Y_mc2[2,1].Intv.hi
@test Y_mc1[1,1].cc == Y_mc2[1,1].cc
@test Y_mc1[1,1].cv == Y_mc2[1,1].cv
@test Y_mc1[2,1].cc == Y_mc2[2,1].cc
@test Y_mc1[2,1].cv == Y_mc2[2,1].cv
@test Y_mc1[1,2].Intv.lo == Y_mc2[1,2].Intv.lo
@test Y_mc1[1,2].Intv.hi == Y_mc2[1,2].Intv.hi
@test Y_mc1[2,2].Intv.lo == Y_mc2[2,2].Intv.lo
@test Y_mc1[2,2].Intv.hi == Y_mc2[2,2].Intv.hi
@test Y_mc1[1,2].cc == Y_mc2[1,2].cc
@test Y_mc1[1,2].cv == Y_mc2[1,2].cv
@test Y_mc1[2,2].cc == Y_mc2[2,2].cc
@test Y_mc1[2,2].cv == Y_mc2[2,2].cv

# Test non-sparse preconditioner
h1 = z_mc2
hj1 = Y_mc2
Y1 = [1.0 2.0; 4.0 5.0; 8.0 9.0]
nx = 2
EAGO.Precondition!(h1,hj1,Y1,nx)
@test h1[1].cc == 7.5
@test h1[1].cv == 1.5
@test h1[1].cc_grad == SVector{2,Float64}([1.5, 0.0])
@test h1[1].cv_grad == SVector{2,Float64}([1.5, 0.0])
@test h1[1].Intv.lo == 0.0
@test h1[1].Intv.hi == 10.5
@test h1[2].cc == 42.5
@test h1[2].cv == 8.5
@test h1[2].cc_grad == SVector{2,Float64}([8.5, 0.0])
@test h1[2].cv_grad == SVector{2,Float64}([8.5, 0.0])
@test h1[2].Intv.lo == 0.0
@test h1[2].Intv.hi == 59.5

@test hj1[1,1].cc == 7.5
@test hj1[1,1].cv == 1.5
@test hj1[1,1].cc_grad == SVector{2,Float64}([1.5, 0.0])
@test hj1[1,1].cv_grad == SVector{2,Float64}([1.5, 0.0])
@test hj1[1,1].Intv.lo == 0.0
@test hj1[1,1].Intv.hi == 10.5
@test hj1[2,1].cc == 42.5
@test hj1[2,1].cv == 8.5
@test hj1[2,1].cc_grad == SVector{2,Float64}([8.5, 0.0])
@test hj1[2,1].cv_grad == SVector{2,Float64}([8.5, 0.0])
@test hj1[2,1].Intv.lo == 0.0
@test hj1[2,1].Intv.hi == 59.5

@test hj1[1,2].cc == 7.5
@test hj1[1,2].cv == 1.5
@test hj1[1,2].cc_grad == SVector{2,Float64}([1.5, 0.0])
@test hj1[1,2].cv_grad == SVector{2,Float64}([1.5, 0.0])
@test hj1[1,2].Intv.lo == 0.0
@test hj1[1,2].Intv.hi == 10.5
@test hj1[2,2].cc == 42.5
@test hj1[2,2].cv == 8.5
@test hj1[2,2].cc_grad == SVector{2,Float64}([8.5, 0.0])
@test hj1[2,2].cv_grad == SVector{2,Float64}([8.5, 0.0])
@test hj1[2,2].Intv.lo == 0.0
@test hj1[2,2].Intv.hi == 59.5

end
