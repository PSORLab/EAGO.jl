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

  EAGO.ExtendDivStorage!(S1,B1,tL,tU,X1,X2,pos,LBDorig,UBDorig)
  EAGO.ExtendDivStorage!(S2,B2,tL,tU,X1,X2,pos,LBDorig,UBDorig)

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

@testset "Ipopt Callback Functions" begin

  cb = EAGO.callback_storage()

  opts = EAGO.EAGO_Inner_NLP(EAGO_NLPSolver())
  opts.f = x -> x[1]*x[4]*(x[1] + x[2] + x[3]) + x[3]
  opts.g = x ->  [x[1]*x[2]*x[3]*x[4]; x[1]^2 + x[2]^2 + x[3]^2 + x[4]^2]
  opts.numVar = 4
  opts.numConstr = 2

  obj_factor = 0.0
  lambda = zeros(Float64,4)

  x = [1.0,2.0,3.0,4.0]
  X = [Interval(0.0,2.0),Interval(1.0,3.0),Interval(2.0,4.0),Interval(3.0,5.0)]
  f_grad = zeros(Float64,4)
  g = zeros(Float64,2)

  rows = zeros(Int32,8)
  cols = zeros(Int32,8)
  values = zeros(Float64,8)
  for i = 1:2
    cols[(4*(i-1)+1):(4*i)] = 1:4
    rows[(4*(i-1)+1):(4*i)] = ones(4)*i
  end
  cb.row_temp_Ipopt_LBD = rows
  cb.col_temp_Ipopt_LBD = cols

  cb.IPOPT_LBD_eval_f = (x::Vector{Float64}, X) -> EAGO.IPOPT_LBD_eval_f(x, X, opts)
  cb.IPOPT_LBD_eval_g = (x::Vector{Float64}, X::Vector{Interval{Float64}}) -> EAGO.IPOPT_LBD_eval_g(x, X, opts)

  fval = EAGO.IPOPT_LBD_eval_f(x,X,opts)
  @test fval == 12.0

  EAGO.IPOPT_LBD_eval_grad_f!(x,X,f_grad,opts)
  EAGO.IPOPT_LBD_eval_g!(x,X,g,opts)
  f_grad1 = copy(f_grad)
  g1 = copy(g)
  @test f_grad1 == Float64[9.00 0.00 1.00 0.00]
  @test g1 == Float64[6.00 36.00]

  gval = EAGO.IPOPT_LBD_eval_g(x,X,opts)
  @test gval == loat64[6.00 36.00]

  mode = :Values
  EAGO.IPOPT_LBD_eval_jac_g!(x,X,mode,rows,cols,values,opts,cb)
  @test values == [6.00, 0.00, 0.00, 0.00, 2.00, 4.00, 6.00, 8.00]

  mode = :Structure
  EAGO.IPOPT_LBD_eval_jac_g!(x,X,mode,rows,cols,values,opts,cb)
  @test values == [6.00, 0.00, 0.00, 0.00, 2.00, 4.00, 6.00, 8.00]
end

end
