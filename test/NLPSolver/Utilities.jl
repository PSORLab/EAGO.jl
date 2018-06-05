module Utilities_NLP_Test

using Compat
using Compat.Test
using EAGO
using JuMP
using MathProgBase

@testset "Repeat Node" begin

  struct testsolver
    max_reduce_rept::Int
    tol_reduce_rept
  end

  struct testopt
    numVar::Int
    solver::testsolver
  end

  opt = testopt(2,testsolver(2,0.3))

  rep = false
  X0a = [Interval(1.0,4.0),Interval(1.0,10.0)]
  X1a = [Interval(1.0,4.0),Interval(1.0,2.0)]
  X2a = copy(X0a)
  @test repeat_DR(X1a,X0a,opt,3,rep) == false
  @test repeat_DR(X1a,X0a,opt,1,rep) == true
  @test repeat_DR(X2a,X0a,opt,1,rep) == false

end

@testset "Extended Division Store" begin

struct testeddiv
  ed_flag::Bool
  box2
end

function ExtendDivStorage!(S::BnBSolver,B::BnBModel{T},tL::Float64,tU::Float64,
                           X1::Vector{V},
                           X2::Vector{V},pos::Int64,
                           LBDorig::Float64,UBDorig::Float64) where {T<:AbstractFloat,V<:AbstractInterval}
    if (S.opt[3].ed_flag)
        push!(B.box,X1,S.opt[3].box2)
        push!(B.LBD,tL,LBDorig)
        push!(B.UBD,tU,UBDorig)
        push!(B.id,B.max_id+1,B.max_id+2)
        push!(B.pos,pos+1,pos+1)
        B.max_id += 2
    else
        push!(B.box,X1,X2)
        push!(B.LBD,tL,tL)
        push!(B.UBD,tU,tU)
        push!(B.id,B.max_id+1,B.max_id+2)
        push!(B.pos,pos+1,pos+1)
        B.max_id += 2
    end
end

  S1 = BnBSolver()
  S2 = BnBSolver()
  S1.opt = [1,2,testeddiv(true,[Interval(-3.0,-1.0),Interval(-8.0,-2.0)])]
  S2.opt = [1,2,testeddiv(false)]
  B1 = BnBModel([Interval(1.0,3.0),Interval(2.0,8.0)])
  B2 = BnBModel([Interval(1.0,3.0),Interval(2.0,8.0)])
  tL = -1.0
  tU = 2.0
  X1 = [Interval(1.0,3.0),Interval(2.0,8.0)]
  X2 = [Interval(1.0,5.0),Interval(3.0,7.0)]
  pos = 1
  LBDorig = -2.0
  UBDorig = 4.0

  ExtendDivStorage!(S1,B1,tL,tU,X1,X2,pos,LBDorig,UBDorig)
  ExtendDivStorage!(S2,B2,tL,tU,X1,X2,pos,LBDorig,UBDorig)

  @test B1.box[end-1] == X1
  @test B1.box[end] == [Interval(-3.0,-1.0),Interval(-8.0,-2.0)]
  @test B1.LBD[end] == LBDorig
  @test B1.LBD[end-1] == tL
  @test B1.UBD[end] == UBDorig
  @test B1.UBD[end-1] == tU
  @test B1.id[end] == 3
  @test B1.id[end-1] == 2
  @test B1.pos[end] == 3
  @test B1.pos[end-1] == 2
  @test B1.max_id == 3

  @test B2.box[end-1] == X1
  @test B2.box[end] == X2
  @test B2.LBD[end-1] == tL
  @test B2.LBD[end] == tL
  @test B2.UBD[end-1] == tU
  @test B2.UBD[end] == tU
  @test B2.id[end] == 3
  @test B2.id[end-1] == 2
  @test B2.pos[end] == 3
  @test B2.pos[end-1] == 2
  @test B2.max_id == 3

end

end
