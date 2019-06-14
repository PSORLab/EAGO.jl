module Check_Node_Access_Functions

    using Compat
    using Compat.Test
    using EAGO

    @testset "Node Access Functions" begin
        x = EAGO.NodeBB(Float64[1.0,5.0], Float64[2.0,6.0], -3.4, 2.1, 2, 1, true)

        @test EAGO.lower_variable_bounds(x) == Float64[1.0,5.0]
        @test EAGO.upper_variable_bounds(x) == Float64[2.0,6.0]
        @test EAGO.lower_bound(x) == -3.4
        @test EAGO.upper_bound(x) == 2.1
        @test EAGO.depth(x) == 2
        @test EAGO.last_branch(x) == 1
        @test EAGO.branch_direction(x) == true
    end

end
