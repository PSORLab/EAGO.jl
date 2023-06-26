@testset "Duality-Based Bound Tightening" begin

    ylower = Float64[1.0, 1.0, 1.0, 1.0]
    yupper = Float64[4.0, 4.0, 4.0, 4.0]
    ymult_lo = Float64[50.0, 0.0, 1.0, 0.0]
    ymult_hi = Float64[0.0, 0.0, 0.8, 3.0]
    isint = Bool[false, false]
    n = EAGO.NodeBB(ylower, yupper, isint, true, -Inf, Inf, 2, 1, 1, EAGO.BD_NONE, 1, 0.1)
    @inferred EAGO.variable_dbbt!(n, ymult_lo, ymult_hi, 1.0, 3.0, 4)
    lvb = n.lower_variable_bounds
    uvb = n.upper_variable_bounds

    @test lvb[1] == 1.0; @test uvb[1] == 1.04
    @test lvb[2] == 1.0; @test uvb[2] == 4.0
    @test lvb[3] == 1.0; @test uvb[3] == 3.0
    # lvb[3] doesn't tighten since dbbt assumes variables aren't fixed
    @test isapprox(lvb[4], 3.33333333, atol= 1E-5); @test uvb[4] == 4.0
end