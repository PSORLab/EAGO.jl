
intype(::Type{IntervalArithmetic.Interval{T}}) where T<:AbstractFloat = T

mutable struct McCormickParamters
  env_max_int::Int64 # Adjusted
  env_tol::Float64 # Adjusted
  mu::Int64 # Adjusted
  valid_check::Bool # Adjusted
  valid_tol::Float64
  multivar_refine::Bool # Adjusted
  mv_tol::Float64 # Adjusted
  outer_rnding::Bool # Adjusted
  outer_param::Float64 # Adjusted
  McCormickParamters() = new(100,
                             1E-10,
                             0,
                             false,
                             1E-8,
                             false,
                             1E-15,
                             false,
                             0.0)
end

const MC_param = McCormickParamters()

flt_un = Union{Float16,Float32,Float64}

#abstract type AbstractMC <: Real end
"""
    SMCg{N,V,T<:AbstractFloat}

`SMCg` is the smooth McCormick (w/ gradient) structure which is used to overload
standard calculations. The fields are:
* `cc::T`: Concave relaxation
* `cv::T`: Convex relaxation
* `cc_grad::SVector{N,T}`: (Sub)gradient of concave relaxation
* `cv_grad::SVector{N,T}`: (Sub)gradient of convex relaxation
* `Intv::V`: Interval bounds
* `cnst::Bool`: Flag for whether the bounds are constant
* `IntvBox::Vector{V}`: Decision space constraints for the affine interval bound tightening
* `xref::Vector{T}`: Reference point for affine interval bound tightening
"""
struct SMCg{N,V<:AbstractInterval,T<:AbstractFloat} <: Real
  cc::T
  cv::T
  cc_grad::SVector{N,T}
  cv_grad::SVector{N,T}
  Intv::V
  cnst::Bool

  function SMCg{N,V,T}(cc1::T,cv1::T,cc_grad1::SVector{N,T},cv_grad1::SVector{N,T},
                Intv1::V,cnst1::Bool) where {N,V<:AbstractInterval,T<:AbstractFloat}
    new(cc1,cv1,cc_grad1,cv_grad1,Intv1,cnst1)
  end
end

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

"""SMC(y::Interval) initializes the differentiable McCormick object with an interval
"""
SMCg{N,V,T}(y::V) where {N,V,T} = SMCg(y.hi,y.lo,[],[],y,true)
SMCg{N,V,T}(val,Intv::V) where {N,V,T} = SMCg(val,val,[],[],Intv,true)
