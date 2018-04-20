"""
    seed_g(T::Type,j::Int64,N::Int64)

Creates a `x::SVector{N,T}` object that is one at `x[j]` and zero everywhere else.
"""
function seed_g(T::Type,j::Q,N::Q) where {Q<:Integer}
    return SVector{N,T}([i == j ? 1.0 : 0.0 for i=1:N])
end

"""
    grad(x::SMCg{N,T},j,n) where {N,T}

sets convex and concave (sub)gradients of length `n` of `x` to be `1` at index `j`
"""
function grad(x::SMCg{N,V,T},j::Q) where {N,V,Q<:Integer,T<:AbstractFloat}
  sv_grad::SVector{N,T} = seed_g(T,j,N)
  return SMCg{N,V,T}(x.cc,x.cv,sv_grad,sv_grad,x.Intv,x.cnst,x.IntvBox,x.xref)
end

"""
    zgrad(x::SMCg{N,T},n::Int64) where {N,T}

sets convex and concave (sub)gradients of length `n` to be zero
"""
function zgrad(x::SMCg{N,V,T}) where {N,V,T<:AbstractFloat}
  grad::SVector{N,T} = @SVector zeros(N)
  return SMCg{N,V,T}(x.cc,x.cv,grad,grad,x.Intv,x.cnst,x.IntvBox,x.xref)
end

"""
    mid3(x::T,y::T,z::T)

Calculates the midpoint of three numbers returning the value and the index.
"""
function mid3(x::T,y::T,z::T) where {T<:AbstractFloat}
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
function mid_grad(cc_grad::SVector{N,T}, cv_grad::SVector{N,T}, id::Q) where {N,Q<:Integer,T<:AbstractFloat}
  if (id == 1)
    return cc_grad
  elseif (id == 2)
    return cv_grad
  elseif (id == 3)
    return zeros(SVector{N,T})
  else
    error("Invalid mid3 position")
  end
end

"""
    line_seg(x0::T,x1::T,y1::T,x2::T,y2::T)

Calculates the value of the line segment between `(x1,y1)` and `(x2,y2)` at `x = x0`.
"""
function line_seg(x0::T,x1::T,y1::T,x2::T,y2::T) where {T<:AbstractFloat}
   if (x2-x1) == zero(T)
     return y1
   else
     return y1*((x2-x0)/(x2-x1)) + y2*((x0-x1)/(x2-x1))
   end
end

"""
    dline_seg(x0::T,x1::T,y1::T,x2::T,y2::T)

Calculates the value of the slope line segment between `(x1,y1)` and `(x2,y2)`
defaults to evaluating the derivative of the function if the interval is tight.
"""
function dline_seg(x0::T,x1::T,y1::T,x2::T,y2::T,d::T) where {T<:AbstractFloat}
    if (x2 == x1)
      return d
    else
      return (y2-y1)/(x2-x1)
    end
end

"""
    grad_calc(cv::T,cc::T,int1::Int64,int2::Int64,dcv::SVector{N,T},dcc::SVector{N,T}) where {N,T}

(Sub)gradient calculation function. Takes the convex gradient, 'cv', the
concave gradient, 'cc', the mid index values 'int1,int2', and the derivative of
the convex and concave envelope functions 'dcv,dcc'.
"""
function grad_calc(cv::SVector{N,T},cc::SVector{N,T},int1::Q,int2::Q,dcv::T,dcc::T) where {N,Q<:Integer,T<:AbstractFloat}
  cv_grad::SVector{N,T} = dcv*( int1==1 ? cv :( int1==2 ? cv : zeros(SVector{N,T})))
  cc_grad::SVector{N,T} = dcc*( int2==1 ? cc :( int2==2 ? cc : zeros(SVector{N,T})))
  return cv_grad, cc_grad
end

"""
    tighten_subgrad(cc,cv,cc_grad,cv_grad,Xintv,Xbox,xref)

Tightens the interval bounds using subgradients. Inputs:
* `cc::T`: concave bound
* `cv::T`: convex bound
* `cc_grad::SVector{N,T}`: subgradient/gradient of concave bound
* `cv_grad::SVector{N,T}`: subgradient/gradient of convex bound
* `Xintv::Interval{T}`: Interval domain of function
* `Xbox::Vector{Interval{T}}`: Original decision variable bounds
* `xref::Vector{T}`: Reference point in Xbox
"""
function tighten_subgrad(cc::T,cv::T,cc_grad::SVector{N,T},cv_grad::SVector{N,T},
                         Xintv::V,Xbox::Vector{V},xref::Vector{T}) where {N,V,T<:AbstractFloat}
  if (length(Xbox)>0 && Xbox[1]!=∅)
    upper_refine::V = convert(V,cc)
    lower_refine::V = convert(V,cv)
    for i=1:N
      upper_refine = upper_refine + cc_grad[i]*(Xbox[i]-xref[i])
      lower_refine = lower_refine + cv_grad[i]*(Xbox[i]-xref[i])
    end
    return max(lower_refine.lo,Xintv.lo), min(upper_refine.hi,Xintv.hi)
  end
end

"""
    outer_rnd!(Intv::Interval{T})

Outer rounds the interval `Intv` by `MC_param.outer_param`.
"""
function outer_rnd(Intv::V) where {V}
    return Intv.lo-MC_param.outer_param, Intv.hi+MC_param.outer_param
end

"""
    isequal(x::S,y::S,atol::T,rtol::T)

Checks that `x` and `y` are equal to with absolute tolerance `atol` and relative
tolerance `rtol`.
"""
function isequal(x::S,y::S,atol::T,rtol::T) where {S,T<:AbstractFloat}
  return (abs(x-y) < (atol + 0.5*abs(x+y)*rtol))
end

"""
    Intv(x::SMCg{N,T})
"""
function Intv(x::SMCg{N,V,T}) where {N,V,T<:AbstractFloat}
  return x.Intv
end

"""
  cut!(x::SMCg{N,V,T})
"""
function cut(xL::T,xU::T,cv::T,cc::T,cv_grad::SVector{N,T},cc_grad::SVector{N,T}) where {N,T<:AbstractFloat}
    if (cc > xU)
      cco::T = xU
      cc_grado::SVector{N,T} = zeros(SVector{N,T})
    else
      cco = cc
      cc_grado = cc_grad
    end
    if (cv < xL)
      cvo::T  = xL
      cv_grado::SVector{N,T} = zeros(SVector{N,T})
    else
      cvo = cv
      cv_grado = cv_grad
  end
  return cvo,cco,cv_grado,cc_grado
end

const eq_tol_flt64 = Float64(1E-12)
const eq_tol_flt32 = Float32(1E-12)
const eq_tol_flt16 = Float16(1E-12)
eq_tol(::Type{Float64}) = eq_tol_flt64
eq_tol(::Type{Float32}) = eq_tol_flt32
eq_tol(::Type{Float16}) = eq_tol_flt16

function neq(x::T,y::T) where {T<:AbstractFloat}
  abs(x-y)>eq_tol(T)
end

const half_flt64 = Float64(0.5)
const half_flt32 = Float32(0.5)
const half_flt16 = Float16(0.5)
half(::Type{Float64}) = half_flt64
half(::Type{Float32}) = half_flt32
half(::Type{Float16}) = half_flt16

const two_flt64 = Float64(2.0)
const two_flt32 = Float32(2.0)
const two_flt16 = Float16(2.0)
two(::Type{Float64}) = two_flt64
two(::Type{Float32}) = two_flt32
two(::Type{Float16}) = two_flt16

const three_flt64 = Float64(3.0)
const three_flt32 = Float32(3.0)
const three_flt16 = Float16(3.0)
three(::Type{Float64}) = three_flt64
three(::Type{Float32}) = three_flt32
three(::Type{Float16}) = three_flt16

# defines square operator
sqr(x::T) where T = x*x
