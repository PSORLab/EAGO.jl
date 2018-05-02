module IntervalTrignometric_Check

using Compat
using Compat.Test
using EAGO

@testset "trig utilities" begin

@test EAGO.two_pi(Float32) == EAGO.pi232
@test EAGO.two_pi(Float16) == EAGO.pi216
@test EAGO.half_pi(Float64) == EAGO.pi_MCinterval(Float64) / 2.0

end

@testset "sin" begin
    @test sin(MCInterval(0.5)) == MCInterval(0.479425538604203, 0.479425538604203)
    @test sin(MCInterval(0.5, 1.67)) == MCInterval(0.479425538604203, 1.0)
    @test sin(MCInterval(1.67, 3.2)) == MCInterval(-0.058374143427580086, 0.9950833498101802)
    @test sin(MCInterval(2.1, 5.6)) == MCInterval(-1.0, 0.8632093666488737)
    @test sin(MCInterval(0.5, 8.5)) == MCInterval(-1.0, 1.0)
    @test sin(MCInterval(0.0, 6.0)) == MCInterval(-1.0, 1.0)
    @test sin(MCInterval(6.0, 7.0)) == MCInterval(-0.27941549819892586, 0.6569865987187891)
end

@testset "cos" begin
    @test cos(MCInterval(0.5)) == MCInterval(0.8775825618903728, 0.8775825618903728)
    @test cos(MCInterval(0.5, 1.67)) == MCInterval(-0.09904103659872801, 0.8775825618903728)
    @test cos(MCInterval(2.1, 5.6)) == MCInterval(-1.0, 0.7755658785102496)
    @test cos(MCInterval(0.5, 8.5)) == MCInterval(-1.0, 1.0)
    @test cos(MCInterval(1.67, 3.2)) == MCInterval(-1.0, -0.09904103659872801)
    @test cos(MCInterval(6.0, 7.0)) == MCInterval(0.7539022543433046, 1.0)
    @test cos(MCInterval(4.0, 6.0)) == MCInterval(-0.6536436208636119, 0.960170286650366)
    @test cos(MCInterval(2.0, 7.0)) == MCInterval(-1.0, 1.0)
end

@testset "tan" begin
    @test tan(MCInterval(0.5,0.5)) == MCInterval(0.5463024898437905, 0.5463024898437905)
    @test tan(MCInterval(-1.5, -1.4)) == MCInterval(-14.101419947171719, -5.797883715482887)
    @test tan(MCInterval(-0.2, -0.1)) == MCInterval(-0.2027100355086725, -0.10033467208545055)
    @test tan(MCInterval(0.1, 0.2)) == MCInterval(0.10033467208545055, 0.2027100355086725)
    @test tan(MCInterval(1.4, 1.5)) == MCInterval(5.797883715482887, 14.101419947171719)
    @test tan(MCInterval(0.5, 9.97)) == entireMCinterval(Float64)                                # Significant failure (REVISIT)
    @test tan(MCInterval(1.67, 3.2)) == MCInterval(-10.047182299210306, 0.058473854459578645)
    #@test tan(MCInterval(0.0, 3.0)) == entireMCinterval(Float64)      # Significant failure (REVISIT)
end


@testset "Inverse trig" begin
    @test asin(MCInterval(1,1)) == MCInterval(pi/2)#pi_interval(Float64)/2
    @test asin(MCInterval(0.9, 2.0)) == asin(MCInterval(0.9, 1.0))
    @test asin(MCInterval(3, 4)) == emptyMCinterval(Float64)

    @test acos(MCInterval(1,1)) == MCInterval(0., 0.)
    @test acos(MCInterval(-2.0, -0.9)) == acos(MCInterval(-1.0, -0.9))
    @test acos(MCInterval(3, 4)) == emptyMCinterval(Float64)

    @test atan(MCInterval(-1,1)) == MCInterval(-pi/4, pi/4)
    @test atan(MCInterval(0,0)) == MCInterval(0.0, 0.0)
end

@testset "Trig" begin

    @test sin(MCInterval(-pi/2, 3pi/2)) == MCInterval(-1, 1)
    @test sin(MCInterval(-100.0, 100.0)) == MCInterval(-1, 1)
    @test cos(MCInterval(-pi/2, 3pi/2)) == MCInterval(-1, 1)
    @test cos(MCInterval(-100.0, 100.0)) == MCInterval(-1, 1)
end

end
