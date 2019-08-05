#abstract type AbstractMC <: Real end
"""
    MC

`MC` is the smooth McCormick (w/ gradient) structure which is used to overload
standard calculations. The fields are:
* `cc::Float64`: Concave relaxation
* `cv::Float64`: Convex relaxation
* `cc_grad::SVector{N,Float64}`: (Sub)gradient of concave relaxation
* `cv_grad::SVector{N,Float64}`: (Sub)gradient of convex relaxation
* `Intv::IntervalType`: Interval bounds
* `cnst::Bool`: Flag for whether the bounds are constant
"""
struct MC{N} <: Real
  cv::Float64
  cc::Float64
  Intv::IntervalType
  cv_grad::SVector{N,Float64}
  cc_grad::SVector{N,Float64}
  cnst::Bool

  function MC{N}(cv1::Float64,cc1::Float64,Intv1::IntervalType,
                 cv_grad1::SVector{N,Float64},cc_grad1::SVector{N,Float64},cnst1::Bool) where {N}
    new(cv1,cc1,Intv1,cv_grad1,cc_grad1,cnst1)
  end
end

"""
    MC(y::IntervalType)

Constructs McCormick relaxation with convex relaxation equal to `y.lo` and
concave relaxation equal to `y.hi`.
"""
MC{N}(y::IntervalType) where N = MC{N}(y.lo,y.hi,y,SVector{N,Float64}(zeros(Float64,N)),SVector{N,Float64}(zeros(Float64,N)),true)
MC{N}(y::Float64) where N = MC{N}(IntervalType(y))
MC{N}(y::T) where {N,T<:AbstractIrrational} = MC{N}(IntervalType(y))
MC{N}(cv::Float64, cc::Float64) where N = MC{N}(cv,cc,IntervalType(cv,cc),SVector{N,Float64}(zeros(Float64,N)),SVector{N,Float64}(zeros(Float64,N)),true)

"""
    MC{N}(val,intv::IntervalType,i::Int)

Constructs McCormick relaxation with convex and concave relaxation equal to
`val`, interval bounds equal to `intv`, and subgradients are the i^th standard
unit vector.
"""
MC{N}(val,Intv::IntervalType,i::Int) where N = MC{N}(val,val,Intv,seed_gradient(Float64,i,N),seed_gradient(Float64,i,N),false)

MC{N}(x::MC{N}) where N = MC{N}(x.cv,x.cc,x.Intv,x.cv_grad,x.cc_grad,x.cnst)
