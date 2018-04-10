
intype(::Type{IntervalArithmetic.Interval{T}}) where T<:AbstractFloat = T

mutable struct McCormickParamters
  env_max_int::Int64 # Adjusted
  env_tol::Float64 # Adjusted
  mu::Int64 # Adjusted
  valid_check::Bool # Adjusted
  valid_tol::Float64
  subgrad_refine ::Bool# Adjusted
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
                             false,
                             1E-15,
                             false,
                             0.0)
end

const MC_param = McCormickParamters()

flt_un = Union{Float16,Float32,Float64}

#abstract type AbstractMC <: Real end
"""
    SMCg{N,T<:AbstractFloat}

`SMCg` is the smooth McCormick (w/ gradient) structure which is used to overload
standard calculations. The fields are:
* `cc::T`: Concave relaxation
* `cv::T`: Convex relaxation
* `cc_grad::SVector{N,T}`: (Sub)gradient of concave relaxation
* `cv_grad::SVector{N,T}`: (Sub)gradient of convex relaxation
* `Intv::Interval{T}`: Interval bounds
* `cnst::Bool`: Flag for whether the bounds are constant
* `IntvBox::Vector{Interval{T}}`: Decision space constraints for the affine interval bound tightening
* `xref::Vector{T}`: Reference point for affine interval bound tightening
Affine interval bound tightening is included in the constructor. As well, as a
validity check options and the cut operator is included if the differentiability
is set to nonsmooth.
"""
struct SMCg{N,V,T<:AbstractFloat} <: Real
  cc::T
  cv::T
  cc_grad::SVector{N,T}
  cv_grad::SVector{N,T}
  Intv::V
  cnst::Bool
  IntvBox::SVector{N,V}
  xref::SVector{N,T}

  function SMCg{N,V,T}(cc1::T,cv1::T,cc_grad1::SVector{N,T},cv_grad1::SVector{N,T},
                Intv1::V,cnst1::Bool,Intv1Box::SVector{N,V},
                xref1::SVector{N,T}) where {N,V,T<:AbstractFloat}
    #=
    can uncomment and enable for debugging.... really onlu useful then.....
    if MC_param.valid_check
      if ((cc1+MC_param.valid_tol)<cv1)
        error("cc must be greater than or equal to cv. cc is $cc1. cv is $cv1")
      elseif ((cc1-MC_param.valid_tol)>Intv1.hi)
        error("cc must be less than or equal to upper interval bound. cc is $cc1. Intv.hi is $(Intv1.hi)")
      elseif ((cv1+MC_param.valid_tol)<Intv1.lo)
        error("cv must be greater than or equal to lower interval bound. cv is $cv1. cv is $(Intv1.lo)")
      end
    end
    =#
    new(cc1,cv1,cc_grad1,cv_grad1,Intv1,cnst1,Intv1Box,xref1)
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
SMCg{N,V,T}(y::V,IntvBox,xref1) where {N,V,T} = SMCg(y.hi,y.lo,[],[],y,true,IntvBox,xref1)
SMCg{N,V,T}(val,Intv::V,IntvBox,xref1) where {N,V,T} = SMCg(val,val,[],[],Intv,true,IntvBox,xref1)
