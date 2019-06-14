module ParamIntvPrecond_Tests

using Compat
using Compat.Test
using IntervalArithmetic
using EAGO, LinearAlgebra

@testset "Preconditioners" begin
    opt1 = parametric_interval_params(:None,:Newton,1E-30,1E-6,Int(2),Int(2),100)
    function h(out,x,xold,p,t)
         out[1] = p[1]
     end
    function hj(out,x,xold,p,t)
        out[1,1] = p[1]
    end
    X = [IntervalType(1.0,2.0),IntervalType(1.0,2.0)]
    P = [IntervalType(1.0,2.0),IntervalType(1.0,3.0)]

    H = fill(Interval(0.0),(1,))
    J = fill(Interval(0.0),(1,1))
    Xold = X
    Y = fill(0.0,(1,1))
    nx = 1
    EAGO.dense_precondition(h, hj, H, J, X, Xold, P, Float64[], Y, nx)
    @test isapprox(H[1].lo, 0.66666, atol=1E-3)
    @test isapprox(H[1].hi, 1.33334, atol=1E-3)
end

end
