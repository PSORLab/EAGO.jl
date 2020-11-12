@testset "Explicit SIP Routine" begin
    # Define semi-infinite program
    f(x) = (1/3)*x[1]^2 + x[2]^2 + x[1]/2
    gSIP(x,p) = (1.0 - (x[1]^2)*(p[1]^2))^2 - x[1]*p[1]^2 - x[2]^2 + x[2]

    x_l = [-1000.0, -1000.0]
    x_u = [1000.0, 1000.0]
    p_l = [0.0]
    p_u = [1.0]

    # Create optimizer for use in solving SIP
    opt = optimizer_with_attributes(EAGO.Optimizer, "cut_max_iterations" => 0,
                                                    "subgrad_tighten" => false,
                                                    "absolute_tolerance" => 1E-4,
                                                    "relative_tolerance" => 1E-4,
                                                    "verbosity" => 0)
    m = JuMP.Model(opt)
    sip_result = explicit_sip_solve(x_l, x_u, p_l, p_u, f, [gSIP], sip_absolute_tolerance = 1E-3)

    @test isapprox(sip_result.lower_bound, 0.19446619886176916, atol = 1E-3)
    @test isapprox(sip_result.upper_bound, 0.19500611503848964, atol = 1E-3)
    @test isapprox(sip_result.xsol[1], -0.7500000115038946, atol = 1E-2)
    @test isapprox(sip_result.xsol[2], -0.6184706298867955, atol = 1E-2)
end
