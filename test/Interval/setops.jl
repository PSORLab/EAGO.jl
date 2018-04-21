module IntervalSetOps_Check

using Compat
using Compat.Test
using EAGO

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
end

end
