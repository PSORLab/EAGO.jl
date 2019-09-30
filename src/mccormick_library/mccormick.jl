module McCormick

using DiffRules: diffrule
using StaticArrays: @SVector, SVector, zeros, ones

import Base: +, -, *, /, convert, in, isempty, one, zero, real, eps, max, min,
             abs, inv, exp, exp2, exp10, expm1, log, log2, log10, log1p, acosh, sech,
             csch, coth, acsch, acoth, asech,
             sqrt, sin, cos, tan, min, max, sec, csc, cot, ^, step, sign, intersect,
             promote_rule, asinh, atanh, tanh, atan, asin, cosh, acos,
             sind, cosd, tand, asind, acosd, atand,
             secd, cscd, cotd, asecd, acscd, acotd

import IntervalArithmetic: dist, mid, pow, +, -, *, /, convert, in, isempty,
                           one, zero, real, eps, max, min, abs, exp,
                           expm1, log, log2, log10, log1p, sqrt, ^,
                           sin, cos, tan, min, max, sec, csc, cot, step,sech,
                           csch, coth, acsch, acoth, asech,
                           sign, dist, mid, pow, Interval, interval, sinh, cosh,
                           ∩, IntervalBox, pi_interval, bisect, isdisjoint, length,
                           atan, asin, acos,
                           sind, cosd, tand, asind, acosd, atand,
                           secd, cscd, cotd, asecd, acscd, acotd, half_pi, setrounding

import IntervalContractors: plus_rev, mul_rev, min_rev, max_rev, minus_rev, div_rev, exp_rev,
       exp2_rev, exp10_rev, expm1_rev, log_rev, log2_rev, log10_rev,
       log1p_rev, sin_rev, cos_rev, tan_rev, asin_rev, acos_rev, atan_rev,
       sinh_rev, cosh_rev, tanh_rev, asinh_rev, acosh_rev, atanh_rev,
       abs_rev, sqr_rev, sqrt_rev, power_rev

import Base.MathConstants.golden

# Export forward operators
export MC, cc, cv, Intv, lo, hi,  cc_grad, cv_grad, cnst, +, -, *, /, convert,
       one, zero, dist, real, eps, mid, exp, exp2, exp10, expm1, log, log2,
       log10, log1p, acosh, sqrt, sin, cos, tan, min, max, sec, csc, cot, ^,
       abs, step, sign, pow, in, isempty, intersect, length,
       acos, asin, atan, sinh, cosh, tanh, asinh, atanh, inv, sqr, sech,
       csch, coth, acsch, acoth, asech,
       sind, cosd, tand, asind, acosd, atand,
       sinhd, coshd, tanhd, asinhd, acoshd, atanhd,
       secd, cscd, cotd, asecd, acscd, acotd,
       secdh, cschd, cothd, asechd, acschd, acothd

# Export inplace operators
export plus!, mult!, min!, max!, minus!, div!, exp!, exp2!, exp10!, expm1!,
       log!, log2!, log10!, log1p!, sin!, cos!, tan!, asin!, acos!, atan!,
       sinh!, cosh!, tanh!, asinh!, acosh!, atanh!, abs!, sqr!, sqrt!, pow!

export seed_gradient, RelaxTag, NS, MV, Diff

# Export reverse operators
export plus_rev, mul_rev, min_rev, max_rev, minus_rev, div_rev, exp_rev,
       exp2_rev, exp10_rev, expm1_rev, log_rev, log2_rev, log10_rev,
       log1p_rev, sin_rev, cos_rev, tan_rev, asin_rev, acos_rev, atan_rev,
       sinh_rev, cosh_rev, tanh_rev, asinh_rev, acosh_rev, atanh_rev,
       abs_rev, sqr_rev, sqrt_rev, power_rev

# Export utility operators
#=
export grad, zgrad, ∩, mid3, MC_param, mid_grad, seed_g, line_seg, dline_seg,
       outer_rnd, cut, set_valid_check, set_subgrad_refine, set_multivar_refine,
       set_outer_rnd, tighten_subgrad, set_iterations, set_tolerance,
       default_options, value, mincv, maxcc, promote_rule
=#
export mc_opts, gen_expansion_params, gen_expansion_params!, implicit_relax_h,
       implicit_relax_h!, implicit_relax_f, implicit_relax_fg


function __init__()
      setrounding(Interval, :accurate)
end

abstract type RelaxTag end
struct NS <: RelaxTag end
struct MV <: RelaxTag end
struct Diff <: RelaxTag end

const MC_ENV_MAX_INT = 100
const MC_ENV_TOL = 1E-10
const MC_DIFF_MU = 1
const MC_MV_TOL = 1E-8

const IntervalConstr = interval
const Half64 = Float64(0.5)
const Two64 = Float64(2.0)
const Three64 = Float64(3.0)
const EqualityTolerance = Float64(1E-12)
const DegToRadIntv = pi_interval(Float64)/Interval(180.0)
const one_intv = one(Interval{Float64})
const half_intv = Interval{Float64}(0.5)
const two_intv = Interval{Float64}(2.0)
const log2_intv = log(Interval{Float64}(2.0))
const log10_intv = log(Interval{Float64}(10.0))

########### number list
int_list = [Int8,UInt8,Int16,UInt16,
            Int32,UInt32,Int64,UInt64,Int128,UInt128]
float_list = [Float16,Float32,Float64]

########### differentiable functions unitary functions
CVList = [:cosh,:exp,:exp2,:exp10] ### function is convex
CCList = [:acosh,:log,:log2,:log10,:sqrt] ### function is concave
CCtoCVList = [:asin,:sinh,:atanh,:tan] ### function is concave then convex
CVtoCCList = [:atan,:acos,:tanh,:asinh] ### function is convex then concave
Template_List = union(CVList,CCList,CCtoCVList,CVtoCCList)

########### non differentiable and non-unitary functions
OtherList = [:sin,:cos,:min,:max,:abs,:step, :sign, :inv, :*, :+, :-, :/,
:promote_rule, :convert, :one, :zero, :real, :dist, :eps, :fma, :^]

"""
    seed_gradient(j::Int64, x::Val{N})

Creates a `x::SVector{N,Float64}` object that is one at `x[j]` and zero everywhere else.
"""
function seed_gradient(j::Int64, x::Val{N}) where N
  return SVector{N,Float64}([i == j ? 1.0 : 0.0 for i = 1:N]...)
end

"""
    mid3(x::T,y::T,z::T)

Calculates the midpoint of three numbers returning the value and the index.
"""
function mid3(x::Float64,y::Float64,z::Float64)
  (((x>=y)&&(y>=z))||((z>=y)&&(y>=x))) && (return y,2)
  (((y>=x)&&(x>=z))||((z>=x)&&(x>=y))) && (return x,1)
  return z,3
end

"""
    mid_grad(cc_grad::SVector{N,T}, cv_grad::SVector{N,T}, id::Int64)

Takes the concave relaxation gradient 'cc_grad', the convex relaxation gradient
'cv_grad', and the index of the midpoint returned 'id' and outputs the appropriate
gradient according to McCormick relaxation rules.
"""
function mid_grad(cc_grad::SVector{N,Float64}, cv_grad::SVector{N,Float64}, id::Int64) where {N, T <: RelaxTag}
  if (id == 1)
    return cc_grad
  elseif (id == 2)
    return cv_grad
  elseif (id == 3)
    return zeros(SVector{N,Float64})
  end
end

"""
    dline_seg(f::Function, df::Function, x::Float64, xL::Float64, xU::Float64)

Calculates the value of the slope line segment between `(xL, f(xL))` and `(xU, f(xU))`
defaults to evaluating the derivative of the function if the interval is tight.
"""
@inline function dline_seg(f::Function, df::Function, x::Float64, xL::Float64, xU::Float64)
    delta = xU - xL
    if delta == 0.0
        return f(x), df(x)
    else
        yL = f(xL)
        yU = f(xU)
        return (yL*(xU - x) + yU*(x - xL))/delta, (yU - yL)/delta
    end
end
@inline function dline_seg(f::Function, df::Function, x::Float64, xL::Float64, xU::Float64, n::Int64)
    delta = xU - xL
    if delta == 0.0
        return f(x, n), df(x, n)
    else
        yL = f(xL, n)
        yU = f(xU, n)
        return (yL*(xU - x) + yU*(x - xL))/delta, (yU - yL)/delta
    end
end
@inline function dline_seg(f::Function, df::Function, x::Float64, xL::Float64, xU::Float64, c::Float64)
    delta = xU - xL
    if delta == 0.0
        return f(x, c), df(x, c)
    else
        yL = f(xL, c)
        yU = f(xU, c)
        return (yL*(xU - x) + yU*(x - xL))/delta, (yU - yL)/delta
    end
end

"""
    grad_calc(cv::T,cc::T,int1::Int64,int2::Int64,dcv::SVector{N,T},dcc::SVector{N,T}) where {N,T}

(Sub)gradient calculation function. Takes the convex gradient, 'cv', the
concave gradient, 'cc', the mid index values 'int1,int2', and the derivative of
the convex and concave envelope functions 'dcv,dcc'.
"""
function grad_calc(cv::SVector{N,Float64},cc::SVector{N,Float64},int1::Int,int2::Int,dcv::Float64,dcc::Float64) where {N, T <: RelaxTag}
  cv_grad::SVector{N,Float64} = dcv*( int1==1 ? cv : ( int1==2 ? cv : zeros(SVector{N,Float64})))
  cc_grad::SVector{N,Float64} = dcc*( int2==1 ? cc : ( int2==2 ? cc : zeros(SVector{N,Float64})))
  return cv_grad, cc_grad
end

"""
    outer_rnd!(Intv::Interval{T})

Outer rounds the interval `Intv` by `MC_param.outer_param`.
"""
outer_rnd(Intv::Interval{Float64}) = Intv.lo-MC_param.outer_param, Intv.hi+MC_param.outer_param

"""
    isequal(x::Interval{Float64},y::Interval{Float64},atol::Float64,rtol::Float64)

Checks that `x` and `y` are equal to with absolute tolerance `atol` and relative
tolerance `rtol`.
"""
isequal(x::Interval{Float64},y::Interval{Float64},atol::Float64,rtol::Float64) = (abs(x-y) < (atol + 0.5*abs(x+y)*rtol))

"""
    cut
"""
function cut(xL::Float64,xU::Float64,
             cv::Float64,cc::Float64,
             cv_grad::SVector{N,Float64},cc_grad::SVector{N,Float64}) where {N}
    if (cc > xU)
      cco::Float64 = xU
      cc_grado::SVector{N,Float64} = zeros(SVector{N,Float64})
    else
      cco = cc
      cc_grado = cc_grad
    end
    if (cv < xL)
      cvo::Float64 = xL
      cv_grado::SVector{N,Float64} = zeros(SVector{N,Float64})
    else
      cvo = cv
      cv_grado = cv_grad
  end
  return cvo,cco,cv_grado,cc_grado
end

function neq(x::Float64,y::Float64)
  abs(x-y)>EqualityTolerance
end

lo(x::Interval{Float64}) = x.lo
hi(x::Interval{Float64}) = x.hi

function step(x::Interval{Float64})
      isempty(x) && return emptyinterval(x)
      xmin::Float64 = ((x.lo) < 0.0) ? 0.0 : 1.0
      xmax::Float64 = ((x.hi) >= 0.0) ? 1.0 : 0.0
      return Interval{Float64}(xmin,xmax)
end


#abstract type AbstractMC <: Real end
"""
    MC

`MC` is the smooth McCormick (w/ gradient) structure which is used to overload
standard calculations. The fields are:
* `cc::Float64`: Concave relaxation
* `cv::Float64`: Convex relaxation
* `cc_grad::SVector{N,Float64}`: (Sub)gradient of concave relaxation
* `cv_grad::SVector{N,Float64}`: (Sub)gradient of convex relaxation
* `Intv::Interval{Float64}`: Interval bounds
* `cnst::Bool`: Flag for whether the bounds are constant
"""
struct MC{N, T <: RelaxTag} <: Real
  cv::Float64
  cc::Float64
  Intv::Interval{Float64}
  cv_grad::SVector{N,Float64}
  cc_grad::SVector{N,Float64}
  cnst::Bool

  function MC{N,T}(cv1::Float64, cc1::Float64, Intv1::Interval{Float64},
                 cv_grad1::SVector{N,Float64} ,cc_grad1::SVector{N,Float64},
                 cnst1::Bool) where {N, T <: RelaxTag}
    new(cv1,cc1,Intv1,cv_grad1,cc_grad1,cnst1)
  end
end

"""
    MC(y::Interval{Float64})

Constructs McCormick relaxation with convex relaxation equal to `y.lo` and
concave relaxation equal to `y.hi`.
"""
function MC{N,T}(y::Interval{Float64}) where {N, T <: RelaxTag}
    MC{N,T}(y.lo, y.hi, y, SVector{N,Float64}(zeros(Float64,N)),
                           SVector{N,Float64}(zeros(Float64,N)), true)
end
MC{N,T}(y::Float64) where {N, T <: RelaxTag} = MC{N,T}(Interval{Float64}(y))
function MC{N,T}(y::Y) where {N, T <: RelaxTag, Y <: AbstractIrrational}
    MC{N,T}(Interval{Float64}(y))
end
function MC{N,T}(cv::Float64, cc::Float64) where {N, T <: RelaxTag}
    MC{N,T}(cv, cc, Interval{Float64}(cv,cc),
          SVector{N,Float64}(zeros(Float64,N)),
          SVector{N,Float64}(zeros(Float64,N)), true)
end
function MC{N,T}(val::Float64, Intv::Interval{Float64}, i::Int64) where {N, T <: RelaxTag}
    MC{N,T}(val, val, Intv, seed_gradient(i,N), seed_gradient(i,N), false)
end
function MC{N,T}(x::MC{N,T}) where {N, T <: RelaxTag}
    MC{N,T}(x.cv, x.cc, x.Intv, x.cv_grad, x.cc_grad, x.cnst)
end


"""
    grad(x::MC{N,T},j::Int) where {N, T <: RelaxTag}

sets convex and concave (sub)gradients of length `n` of `x` to be `1` at index `j`
"""
function grad(x::MC{N,T},j::Int) where {N, T <: RelaxTag}
  sv_grad::SVector{N,Float64} = seed_gradient(T,j,N)
  return MC{N,T}(x.cc,x.cv,sv_grad,sv_grad,x.Intv,x.cnst)
end

"""
    zgrad(x::SMCg{N,T},n::Int64) where {N,T}

sets convex and concave (sub)gradients of length `n` to be zero
"""
function zgrad(x::MC{N,T}) where {N, T <: RelaxTag}
  grad::SVector{N,Float64} = zeros(SVector{N,Float64})
  return MC{N,T}(x.cc,x.cv,grad,grad,x.Intv,x.cnst)
end

Intv(x::MC) = x.Intv
lo(x::MC) = x.Intv.lo
hi(x::MC) = x.Intv.hi
cc(x::MC) = x.cc
cv(x::MC) = x.cv
cc_grad(x::MC) = x.cc_grad
cv_grad(x::MC) = x.cv_grad
cnst(x::MC) = x.cnst
length(x::MC) = length(x.cc_grad)

"""
    newton(x0::T,xL::T,xU::T,f::Function,df::Function,envp1::T,envp2::T)

Defines a local 1D newton method to solve for the root of `f` between the bounds
`xL` and `xU` using `x0` as a starting point. The derivative of `f` is `df`. The
inputs `envp1` and `envp2` are the envelope calculation parameters.
"""
function newton(x0::Float64, xL::Float64, xU::Float64, f::Function, df::Function, envp1::Float64, envp2::Float64)
  dfk::Float64 = 0.0

  xk::Float64 = max(xL,min(x0,xU))
  fk::Float64 = f(xk,envp1,envp2)

  for i=1:MC_ENV_MAX_INT
    dfk = df(xk, envp1, envp2)
    if (abs(fk) < MC_ENV_TOL)
      return (xk,false)
    end
    (dfk == 0.0) && return (0.0,true)
    if (xk == xL && fk/dfk > 0.0)
      return (xk,false)
    elseif (xk == xU && fk/dfk < 0.0)
      return (xk,false)
    end
    xk = max(xL,min(xU,xk-fk/dfk))
    fk = f(xk,envp1,envp2)
  end
  (0.0,true)
end


"""
    secant(x0::T,x1::T,xL::T,xU::T,f::Function,envp1::T,envp2::T)  where {T<:Real}

Defines a local 1D secant method to solve for the root of `f` between
the bounds `xL` and `xU` using `x0` and `x1` as a starting points. The inputs
`envp1` and `envp2` are the envelope calculation parameters.
"""
function secant(x0::Float64, x1::Float64, xL::Float64, xU::Float64, f::Function, envp1::Float64, envp2::Float64)
  xkm::Float64 = max(xL,min(xU,x0))
  xk::Float64 = max(xL,min(xU,x1))
  fkm::Float64 = f(xkm,envp1,envp2)

  for i=1:MC_ENV_MAX_INT
    fk = f(xk,envp1,envp2)
    Bk::Float64 = (fk-fkm)/(xk-xkm)
    if (abs(fk) < MC_ENV_TOL)
      return (xk,false)
    end
    (Bk == 0.0) && return (0.0, true)
    if ((xk == xL) && (fk/Bk > 0.0))
      return (xk,false)
    elseif ((xk == xU) && (fk/Bk < 0.0))
      return (xk,false)
    end
    xkm = xk
    fkm = fk
    xk = max(xL,min(xU,xk-fk/Bk))
  end
  (0.0,true)
end


"""
    golden_section(xL::T,xU::T,f::Function,envp1::T,envp2::T) where {T<:Real}

Defines a local 1D golden section method to solve for the root of `f` between
the bounds `xL` and `xU` using `x0` as a starting point. Define iteration used
in golden section method. The inputs `envp1` and `envp2` are the envelope
calculation parameters.
"""
function golden_section(xL::Float64,xU::Float64,f::Function,envp1::Float64,envp2::Float64)
  fL::Float64 = f(xL,envp1,envp2)
  fU::Float64 = f(xU,envp1,envp2)

  (fL*fU > 0.0) && error("GOLDEN EXCEPTION")
  xm::Float64 = xU-(2.0-golden)*(xU-xL)
  fm::Float64 = f(xm,envp1,envp2)
  return golden_section_it(1,xL,fL,xm,fm,xU,fU,f,envp1,envp2)
end
"""
    golden_section_it(init::Int64,a::T,fa::T,b::T,fb::T,c::T,
                      fc::T,f::Function,envp1::T,envp2::T) where {T<:Real}

Define iteration used in golden section method. The inputs `fa`,`fb`, and `fc`,
are the function `f` evaluated at `a`,`b`, and `c` respectively. The inputs
`envp1` and `envp2` are the envelope calculation parameters. The value `init` is
the iteration number of the golden section method.
"""
function golden_section_it(init::Int,a::Float64,fa::Float64,b::Float64,fb::Float64,c::Float64,
                                   fc::Float64,f::Function,envp1::Float64,envp2::Float64)
  b_t_x::Bool = (c-b > b-a)
  if (b_t_x)
    x::Float64 = b + (2.0-golden)*(c-b)
  else
    x = b - (2.0-golden)*(b-a)
  end
  itr::Int = init
  if (abs(c-a)<MC_ENV_TOL*(abs(b)+abs(x)) || (itr>MC_ENV_MAX_INT))
    return (c+a)/2.0
  end
  itr += 1
  fx::Float64 = f(x,envp1,envp2)
  if (b_t_x)
    if (fa*fx < 0.0)
      golden_section_it(itr,a,fa,b,fb,x,fx,f,envp1,envp2)
    else
      golden_section_it(itr,b,fb,x,fx,c,fc,f,envp1,envp2)
    end
  else
    if (fa*fb<(fa))
      golden_section_it(itr,a,fa,x,fx,b,fb,f,envp1,envp2)
    else
      golden_section_it(itr,x,fx,b,fb,c,fc,f,envp1,envp2)
    end
  end
end

include("forward_operators/forward.jl")
include("reverse_operators/reverse.jl")
include("implicit_routines/implicit.jl")

end
