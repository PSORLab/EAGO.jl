#=
struct Hybrid_Options{V,T}
    box::SVector{V}
    ref::SVector{T}
end


const hybrid_opts = Hybrid_Options()


"""
    box

Access box bounds for original problem.
"""
box(x::HybridMC) = x.box

"""
    ref

Access reference point for original problem.
"""
ref(x::HybridMC) = x.reg

"""
    boxlo

Access lower box bounds for original problems.
"""
boxlo(x.HybridMC) = lo(x.box)
"""
    boxhi

Access upper box bounds for original problems.
"""
boxhi(x.HybridMC) = hi(x.box)


###################### Subgradient tightening function #########################
"""
    Tighten_Subgrad
"""
function Tighten_Subgrad(x::HybridSMC{N,V,T}) where {N,V,T}
    xIntv = Intv(x)
    lower_cut::T = cc_grad(x)*(ref()-boxlo())
    upper_cut::T = cv_grad(x)*(ref()-boxhi())
    if (lower_cut > xIntv.lo)
        if (upper_cut < xIntv.hi)
            return HybridSMC{N,V,T}(SMCg{N,V,T}(cc(x),cv(x),cc_grad(x),cv_grad(x),V(lower_cut,upper_cut)))
        else
            return HybridSMC{N,V,T}(SMCg{N,V,T}(cc(x),cv(x),cc_grad(x),cv_grad(x),V(lower_cut,xIntv.hi)))
        end
    else
        if (upper_cut < xIntv.hi)
            return HybridSMC{N,V,T}(SMCg{N,V,T}(cc(x),cv(x),cc_grad(x),cv_grad(x),V(xIntv.lo,upper_cut)))
        else
            return HybridSMC{N,V,T}(SMCg{N,V,T}(cc(x),cv(x),cc_grad(x),cv_grad(x),Intv(x)),box(x),ref(x))
        end
    end
end
=#
