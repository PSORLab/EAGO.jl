module Subgradrefine_Test

using Compat
using Compat.Test
using EAGO
using IntervalArithmetic
using StaticArrays

@testset "HybridMC Type & Utilities" begin
    box = SVector{2,Interval{Float64}}([Interval(-1.0,5.0),Interval(-6.5,8.3)])
    ref = SVector{2,Float64}([1.0,2.0])
    hopt = EAGO.Hybrid_Options{2,Interval{Float64},Float64}(box,ref,true)
    @test EAGO.box(hopt) == box
    @test EAGO.ref(hopt) == ref
    @test hopt.sub_on == true

    set_hybrid_box!(box,ref,true)
    @test EAGO.hybrid_opts[1].box == box
    @test EAGO.hybrid_opts[1].ref == ref
    @test EAGO.hybrid_opts[1].sub_on == true

    @test EAGO.boxlo(hopt) == [box[i].lo for i=1:length(box)]
    @test EAGO.boxhi(hopt) == [box[i].hi for i=1:length(box)]

    N = 2
    V = Interval{Float64}
    T = Float64
    cci = 2.0
    cvi = -1.0
    cc_gradi = SVector{N,T}([1.0,0.0])
    cv_gradi = SVector{N,T}([1.0,0.0])
    Intvi = Interval(-4.0,3.0)

    x = HybridMC{N,V,T}(SMCg{N,V,T}(cci,cvi,cc_gradi,cv_gradi,Intvi,false))
    @test cc(x) == cci
    @test cv(x) == cvi
    @test cc_grad(x) == cc_gradi
    @test cv_grad(x) == cv_gradi
    @test Intv(x) == Intvi
end

@testset "Unitary HybridMC Operators" begin

    N = 2
    V = Interval{Float64}
    T = Float64
    cci = 2.0
    cvi = -1.0
    cc_gradi = SVector{N,T}([1.0,0.0])
    cv_gradi = SVector{N,T}([1.0,0.0])
    Intvi = Interval(1.0,3.0)

    y = HybridMC{N,V,T}(SMCg{N,V,T}(cci,cvi,cc_gradi,cv_gradi,Intvi,false))
    @test one(y) == HybridMC(one(SMCg{N,V,T}(cci,cvi,cc_gradi,cv_gradi,Intvi,false)))
    @test zero(y) == HybridMC(zero(SMCg{N,V,T}(cci,cvi,cc_gradi,cv_gradi,Intvi,false)))
    @test real(y) == HybridMC(real(SMCg{N,V,T}(cci,cvi,cc_gradi,cv_gradi,Intvi,false)))
    @test +(y) == y
    @test -(y) == HybridMC{N,V,T}(-SMCg{N,V,T}(cci,cvi,cc_gradi,cv_gradi,Intvi,false))
    @test +(y,y) == 2.0*y
    @test +(1.5,y) == +(y,1.5)

    @test promote_rule(HybridMC{3,Interval{Float64},Float64},Int64) == HybridMC{3,Interval{Float64},Float64}
    @test promote_rule(HybridMC{3,Interval{Float64},Float64},Float64) == HybridMC{3,Interval{Float64},Float64}
    @test promote_rule(HybridMC{3,Interval{Float64},Float64},Interval{Float64}) == HybridMC{3,Interval{Float64},Float64}

    box = SVector{2,Interval{Float64}}([Interval(-0.1,0.2),Interval(-0.1,0.2)])
    ref = SVector{2,Float64}([2.0,2.0])

    set_hybrid_box!(box,ref,false)
    x = HybridMC{N,V,T}(SMCg{N,V,T}(cci,cvi,cc_gradi,cv_gradi,Intvi,false))

    cci1 = 2.0
    cvi1 = -1.0
    cc_gradi1 = SVector{N,T}([2.0,0.0])
    cv_gradi1 = SVector{N,T}([1.0,0.0])
    Intvi1 = Interval(-1.2,2.4)
    x1 = HybridMC{N,V,T}(SMCg{N,V,T}(cci1,cvi1,cc_gradi,cv_gradi,Intvi1,false))

    for (i,z) in enumerate(EAGO.Hybrid_Ops1)
        if (:acosh == z)
            #@test_nowarn eval(Expr(:call,z,x1))
        elseif (:tan == z || :asin == z || :acos == z ||
                :atanh == z || :acosh == z)
        else
            @test_nowarn eval(Expr(:call,z,x))
        end
    end

    set_hybrid_box!(box,ref,true)
    #=
    for (i,z) in enumerate(EAGO.Hybrid_Ops1)
        if (:acosh == z)
            @test_nowarn eval(Expr(:call,z,x1))
        elseif (:tan == z)
        else
            @test_nowarn eval(Expr(:call,z,x))
        end
    end
    =#
    set_hybrid_box!(box,ref,false)

    box = SVector{2,Interval{Float64}}([Interval(-0.2,0.2),Interval(-0.2,0.2)])
    ref = SVector{2,Float64}([0.0,0.0])
    set_hybrid_box!(box,ref,true)
    @test isapprox(Intv.(Tighten_Subgrad(y)).lo,1.000,atol=1E-4)
    @test isapprox(Intv.(Tighten_Subgrad(y)).hi,2.20001,atol=1E-4)
    set_hybrid_box!(box,ref,false)
end

@testset "HybridMC & Number Operators" begin
end

@testset "Binary HybridMC Operators" begin
end

end
