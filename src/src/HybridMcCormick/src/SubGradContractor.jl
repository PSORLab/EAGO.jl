mutable struct Hybrid_Options{N,V<:AbstractInterval,T<:AbstractFloat}
    box::SVector{N,V}
    ref::SVector{N,T}
    sub_on::Bool
end

const hybrid_opts = Hybrid_Options[]


function set_hybrid_box!(box::SVector{N,V},ref::SVector{N,T},subgrad::Bool) where {N,V,T}
    if (isempty(hybrid_opts))
        push!(hybrid_opts,Hybrid_Options{N,V,T}(box,ref,subgrad))
    else
        hybrid_opts[1] = Hybrid_Options{N,V,T}(box,ref,subgrad)
    end
end


"""
    box

Access box bounds for original problem.
"""
box(x::Hybrid_Options) = x.box

"""
    ref

Access reference point for original problem.
"""
ref(x::Hybrid_Options) = x.ref

"""
    boxlo

Access lower box bounds for original problems.
"""
boxlo(x::Hybrid_Options) = [x.box[i].lo for i=1:length(x.box)]
"""
    boxhi

Access upper box bounds for original problems.
"""
boxhi(x::Hybrid_Options) = [x.box[i].hi for i=1:length(x.box)]


###################### Subgradient tightening function #########################
"""
    Tighten_Subgrad
"""
function Tighten_Subgrad(x::HybridMC{N,V,T}) where {N,V,T}
    upper_refine::V = convert(V,cc)
    lower_refine::V = convert(V,cv)
    for i=1:N
      upper_refine = upper_refine + cc_grad[i]*(hybrid_opts[1].box[i]-hybrid_opts[1].ref[i])
      lower_refine = lower_refine + cv_grad[i]*(hybrid_opts[1].box[i]-hybrid_opts[1].ref[i])
    end
    if (lower_refine.lo > xIntv.lo)
        if (upper_refine.hi < xIntv.hi)
            return HybridMC{N,V,T}(SMCg{N,V,T}(cc(x),cv(x),cc_grad(x),cv_grad(x),V(lower_refine.lo,upper_refine)))
        else
            return HybridMC{N,V,T}(SMCg{N,V,T}(cc(x),cv(x),cc_grad(x),cv_grad(x),V(lower_refine.lo,xIntv.hi)))
        end
    else
        if (upper_refine.hi < xIntv.hi)
            return HybridMC{N,V,T}(SMCg{N,V,T}(cc(x),cv(x),cc_grad(x),cv_grad(x),V(xIntv.lo,upper_refine)))
        else
            return x
        end
    end
end

#=
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
=#
