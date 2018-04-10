workspace()

using IntervalArithmetic
using EAGODomainReduction

Y = Interval{Float64}[Interval(1.0,4.0),Interval(1.0,4.0),Interval(1.0,4.0),Interval(1.0,4.0)]
ymult_lo = [50, 1.0, 2.0, 3.0]
ymult_hi = [0, 1.0, 2.0, 3.0]
yLBD = 1.0
yUBD = 3.0
Variable_DR!(Y,ymult_lo,ymult_hi,yLBD,yUBD)
println("Y post full duality-based bound tightening: ", Y)

X = Interval{Float64}[Interval(1.0,4.0),Interval(1.0,4.0),Interval(1.0,4.0),Interval(1.0,4.0)]
mult_lo = [50, 1.0]
mult_hi = [-1.0, 0]
LBD = 1.0
UBD = 2.0
nx = 2
Variable_DR_Imp!(X,mult_lo,mult_hi,LBD,UBD,nx)
println("X post implicit duality-based bound tightening: ", X)
