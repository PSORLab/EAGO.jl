module IntervalSetOps_Check

using Compat
using Compat.Test
using EAGO
using IntervalArithmetic

@test in(1.0,MCInterval(0.0,1.0))
@test ~in(2.0,MCInterval(0.0,1.0))
@test ~in(-1.0,MCInterval(0.0,1.0))
@test ~in(Inf,MCInterval(0.0,1.0))
@test ~in(-Inf,MCInterval(0.0,1.0))

@testset "setdiff" begin
    x = MCInterval(2.0,4.0)
    y = MCInterval(3.0,5.0)

    d = setdiff(x, y)

    @test typeof(d) == Vector{MCInterval{Float64}}
    @test length(d) == 1
    @test d == [MCInterval(2.0,3.0)]
    @test setdiff(y, x) == [MCInterval(4.0,5.0)]

    x = MCInterval(2.0,4.0)
    y = MCInterval(2.0,5.0)

    @test typeof(d) == Vector{MCInterval{Float64}}
    @test length(setdiff(x, y)) == 0
    @test setdiff(y, x) == [MCInterval(4.0,5.0)]

    x = MCInterval(2.0,5.0)
    y = MCInterval(3.0,4.0)
    @test setdiff(x, y) ==[MCInterval(2.0,3.0), MCInterval(4.0,5.0)]

    @test isinterior(y,x)
    @test isinterior(emptyMCinterval(Float64),x)
end

@testset "conversion utils" begin

v64 = flttoMCI(Float64(1.5))
v32 = flttoMCI(Float32(1.5))
v16 = flttoMCI(Float16(1.5))

@test v64.lo == v64.hi == 1.5
@test v32.lo == v32.hi == 1.5
@test v16.lo == v16.hi == 1.5
@test typeof(v64) == MCInterval{Float64}
@test typeof(v32) == MCInterval{Float32}
@test typeof(v16) == MCInterval{Float16}

end

@testset "Interval Checks" begin
x = MCInterval(1.0,Inf)
y = MCInterval(1.0,2.0)
z = MCInterval(1.0,1.0)

@test isthin(z)
@test ~isthin(y)
@test isfinite(y)
@test ~isfinite(x)
@test EAGO.isunbounded(x)

@test EAGO.entireinterval(Float64,MCInterval{Float64}) == entireMCinterval(Float64)
@test EAGO.entireinterval(Float64,Interval{Float64}) == IntervalArithmetic.entireinterval(Float64)
end

end
