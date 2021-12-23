@testset "SIP Res" begin
    # Define semi-infinite program
    f(x) = (1/3)*x[1]^2 + x[2]^2 + x[1]/2
    gSIP(x,p) = (1.0 - (x[1]^2)*(p[1]^2))^2 - x[1]*p[1]^2 - x[2]^2 + x[2]

    x_l = Float64[-1000.0, -1000.0]
    x_u = Float64[1000.0, 1000.0]
    p_l = Float64[0.0]
    p_u = Float64[1.0]

    sip_result = sip_solve(SIPRes(), x_l, x_u, p_l, p_u, f, Any[gSIP], abs_tolerance = 1E-3)

    @test isapprox(sip_result.lower_bound, 0.19446619886176916, atol = 1E-3)
    @test isapprox(sip_result.upper_bound, 0.19500611503848964, atol = 1E-3)
    @test isapprox(sip_result.xsol[1], -0.7500000115038946, atol = 1E-2)
    @test isapprox(sip_result.xsol[2], -0.6184706298867955, atol = 1E-2)
end
#=
@testset "SIP ResRev" begin
    # Define semi-infinite program
    f(x) = (1/3)*x[1]^2 + x[2]^2 + x[1]/2
    gSIP(x,p) = (1.0 - (x[1]^2)*(p[1]^2))^2 - x[1]*p[1]^2 - x[2]^2 + x[2]

    x_l = Float64[-1000.0, -1000.0]
    x_u = Float64[1000.0, 1000.0]
    p_l = Float64[0.0]
    p_u = Float64[1.0]

    sip_result = sip_solve(SIPResRev(), x_l, x_u, p_l, p_u, f, Any[gSIP], abs_tolerance = 1E-3)

    @test isapprox(sip_result.lower_bound, 0.19446619886176916, atol = 1E-3)
    @test isapprox(sip_result.upper_bound, 0.19500611503848964, atol = 1E-3)
    @test isapprox(sip_result.xsol[1], -0.7500000115038946, atol = 1E-2)
    @test isapprox(sip_result.xsol[2], -0.6184706298867955, atol = 1E-2)
end
=#
@testset "SIP Hybrid" begin
    # Define semi-infinite program
    f(x) = (1/3)*x[1]^2 + x[2]^2 + x[1]/2
    gSIP(x,p) = (1.0 - (x[1]^2)*(p[1]^2))^2 - x[1]*p[1]^2 - x[2]^2 + x[2]

    x_l = Float64[-1000.0, -1000.0]
    x_u = Float64[1000.0, 1000.0]
    p_l = Float64[0.0]
    p_u = Float64[1.0]

    sip_result = sip_solve(SIPHybrid(), x_l, x_u, p_l, p_u, f, Any[gSIP], abs_tolerance = 1E-3)

    @test isapprox(sip_result.lower_bound, 0.19446619886176916, atol = 1E-3)
    @test isapprox(sip_result.upper_bound, 0.19500611503848964, atol = 1E-3)
    @test isapprox(sip_result.xsol[1], -0.7500000115038946, atol = 1E-2)
    @test isapprox(sip_result.xsol[2], -0.6184706298867955, atol = 1E-2)
end