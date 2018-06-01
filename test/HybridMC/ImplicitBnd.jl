module ImplicitBnd

using Compat
using Compat.Test
using StaticArrays
using IntervalArithmetic
using EAGO


@testset "Test Base Implicit Routines" begin
    EAGO.set_diff_relax(0)
    opts1 = mc_opts{Float64}(0.5,1,:Dense,:Newton,1,1,1E-10)

    f(x,p) = x[1]*p[1]+p[1]
    g(x,p) = [x[1]*p[1]+p[1];
              x[1]*p[1]+2*p[1]]
    function h1(x,p)
        t1 = x[1]^2
        t2 = x[1]*p[1]
        t3 = 4.0
        t4 = t1 + t2
        t5 = t4 + t3
        return [t5]
    end
    hj1(x,p) = [2*x[1]+p[1]]
    P = [Interval(6.0,9.0)]
    X = [Interval(-0.78,-0.4)]
    p = [7.5]
    pmid = mid.(P)

    xIntv1 = Interval(1.0,3.0)
    xIBox = SVector{1,Interval{Float64}}([xIntv1])
    mBox = mid.(xIBox)
    np = 1
    szero = @SVector zeros(np)
    sone = @SVector ones(np)
    p_mc = [HybridMC{np,Interval{Float64},Float64}(SMCg{np,Interval{Float64},Float64}(p[i],p[i],sone,sone,Interval(P[i].lo,P[i].hi),false)) for i=1:np]
    println("ran to 1")
    param = GenExpansionParams(h1,hj1,X,P,pmid,opts1)
    println("ran to 2")
    hbnds = MC_impRelax(h1,hj1,p_mc,pmid,X,P,opts1,param)
    println("ran to 3")
    fbnds = impRelax_f(f,h1,hj1,X,P,p,pmid,opts1,param)
    println("ran to 4")
    fgbnds = impRelax_fg(f,g,h1,hj1,X,P,p,pmid,opts1,param)
    println("ran to 5")
    param = GenExpansionParams(h1,hj1,X,P,pmid,opts1)
    println("ran to 6")
    @test isapprox(param[2][1].SMC.cc,-0.4,atol=1E-4)
    @test isapprox(param[2][1].SMC.cv,-0.78,atol=1E-4)
    hbnds = MC_impRelax(h1,hj1,p_mc,pmid,X,P,opts1,param)
    @test isapprox(hbnds[1].SMC.cc,-0.5209127114919797,atol=1E-4)
    @test isapprox(hbnds[1].SMC.cv,-0.6307841683146562,atol=1E-4)
    fbnds = impRelax_f(f,h1,hj1,X,P,p,pmid,opts1,param)
    @test isapprox(fbnds.SMC.cc,3.774523731048122,atol=1E-4)
    @test isapprox(fbnds.SMC.cv,2.5572882333553046,atol=1E-4)
    fgbnds = impRelax_fg(f,g,h1,hj1,X,P,p,pmid,opts1,param)
    @test isapprox(fgbnds[1].SMC.cc,3.774523731048122,atol=1E-4)
    @test isapprox(fgbnds[1].SMC.cv,2.5572882333553046,atol=1E-4)
    @test isapprox(fgbnds[2][1].SMC.cc,3.774523731048122,atol=1E-4)
    @test isapprox(fgbnds[2][1].SMC.cv,2.5572882333553046,atol=1E-4)
    @test isapprox(fgbnds[2][2].SMC.cc,11.274523731048122,atol=1E-4)
    @test isapprox(fgbnds[2][2].SMC.cv,10.057288233355305,atol=1E-4)
end

end
