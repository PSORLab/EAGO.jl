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
    cc = 2.0
    cv = -1.0
    cc_grad = SVector{N,T}([1.0,0.0])
    cv_grad = SVector{N,T}([1.0,0.0])
    Intv = Interval(-4.0,3.0)
    #=
    x = HybridMC{N,V,T}(SMCg{N,V,T}(cc,cv,cc_grad,cv_grad))
    @test cc(x) == cc
    @test cv(x) == cv
    @test cc_grad(x) == cc_grad
    @test cv_grad(x) == cv_grad
    @test Intv(x) == Intv
=#
end

@testset "Unitary HybridMC Operators" begin
end

@testset "HybridMC & Number Operators" begin
end

@testset "Binary HybridMC Operators" begin
end

end
