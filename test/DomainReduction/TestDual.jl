module TestDual

using Compat
using Compat.Test
using IntervalArithmetic
using EAGO

@testset "Test Full DR" begin
Y = Interval{Float64}[Interval(1.0,4.0),Interval(1.0,4.0),Interval(1.0,4.0),Interval(1.0,4.0)]
ymult_lo = [50, 1.0, 2.0, 3.0]
ymult_hi = [0, 1.0, 2.0, 3.0]
yLBD = 1.0
yUBD = 3.0
Variable_DR!(Y,ymult_lo,ymult_hi,yLBD,yUBD)
@test 3.95999-1E-4 <= Y[1].lo <= 3.95999+1E-4
@test 4.0-1E-4 <= Y[1].hi <= 4.0+1E-4
@test 2.0-1E-4 <= Y[2].lo <= 2.0+1E-4
@test 4.0-1E-4 <= Y[2].hi <= 4.0+1E-4
@test 3.0-1E-4 <= Y[3].lo <= 3.0+1E-4
@test 4.0-1E-4 <= Y[3].hi <= 4.0+1E-4
@test 3.33333-1E-4 <= Y[4].lo <= 3.33333+1E-4
@test 4.0-1E-4 <= Y[4].hi <= 4.0+1E-4

Y = MCInterval{Float64}[MCInterval(1.0,4.0),MCInterval(1.0,4.0),MCInterval(1.0,4.0),MCInterval(1.0,4.0)]
ymult_lo = [50, 1.0, 2.0, 3.0]
ymult_hi = [0, 1.0, 2.0, 3.0]
yLBD = 1.0
yUBD = 3.0
Variable_DR!(Y,ymult_lo,ymult_hi,yLBD,yUBD)
@test 3.95999-1E-4 <= Y[1].lo <= 3.95999+1E-4
@test 4.0-1E-4 <= Y[1].hi <= 4.0+1E-4
@test 2.0-1E-4 <= Y[2].lo <= 2.0+1E-4
@test 4.0-1E-4 <= Y[2].hi <= 4.0+1E-4
@test 3.0-1E-4 <= Y[3].lo <= 3.0+1E-4
@test 4.0-1E-4 <= Y[3].hi <= 4.0+1E-4
@test 3.33333-1E-4 <= Y[4].lo <= 3.33333+1E-4
@test 4.0-1E-4 <= Y[4].hi <= 4.0+1E-4
end

@testset "Test Implicit DR" begin
X = Interval{Float64}[Interval(1.0,4.0),Interval(1.0,4.0),Interval(1.0,4.0),Interval(1.0,4.0)]
mult_lo = [50, 1.0]
mult_hi = [-1.0, 0]
LBD = 1.0
UBD = 2.0
nx = 2
Variable_DR_Imp!(X,mult_lo,mult_hi,LBD,UBD,nx)
@test 1.0-1E-4 <= X[1].lo <= 1.0+1E-4
@test 4.0-1E-4 <= X[1].hi <= 4.0+1E-4
@test 1.0-1E-4 <= X[2].lo <= 1.0+1E-4
@test 4.0-1E-4 <= X[2].hi <= 4.0+1E-4
@test 3.97999-1E-4 <= X[3].lo <= 3.97999+1E-4
@test 4.0-1E-4 <= X[3].hi <= 4.0+1E-4
@test 3.0-1E-4 <= X[4].lo <= 3.0+1E-4
@test 4.0-1E-4 <= X[4].hi <= 4.0+1E-4

X = MCInterval{Float64}[MCInterval(1.0,4.0),MCInterval(1.0,4.0),MCInterval(1.0,4.0),MCInterval(1.0,4.0)]
mult_lo = [50, 1.0]
mult_hi = [-1.0, 0]
LBD = 1.0
UBD = 2.0
nx = 2
Variable_DR_Imp!(X,mult_lo,mult_hi,LBD,UBD,nx)
@test 1.0-1E-4 <= X[1].lo <= 1.0+1E-4
@test 4.0-1E-4 <= X[1].hi <= 4.0+1E-4
@test 1.0-1E-4 <= X[2].lo <= 1.0+1E-4
@test 4.0-1E-4 <= X[2].hi <= 4.0+1E-4
@test 3.97999-1E-4 <= X[3].lo <= 3.97999+1E-4
@test 4.0-1E-4 <= X[3].hi <= 4.0+1E-4
@test 3.0-1E-4 <= X[4].lo <= 3.0+1E-4
@test 4.0-1E-4 <= X[4].hi <= 4.0+1E-4
end

end
