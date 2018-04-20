"""
"""
function Tighten_Subgrad(x::HybridSMC{N,V,T}) where {N,V,T}
    xIntv = Intv(x)
    lower_cut::T = cc_grad(x)*(ref(x)-boxlo(x))
    upper_cut::T = cv_grad(x)*(ref(x)-boxhi(x))
    if (lower_cut > xIntv.lo)
        if (upper_cut < xIntv.hi)
            return HybridSMC{N,V,T}(SMCg{N,V,T}(cc(x),cv(x),cc_grad(x),cv_grad(x),V(lower_cut,upper_cut)),box(x),ref(x))
        else
            return HybridSMC{N,V,T}(SMCg{N,V,T}(cc(x),cv(x),cc_grad(x),cv_grad(x),V(lower_cut,xIntv.hi)),box(x),ref(x))
        end
    else
        if (upper_cut < xIntv.hi)
            return HybridSMC{N,V,T}(SMCg{N,V,T}(cc(x),cv(x),cc_grad(x),cv_grad(x),V(xIntv.lo,upper_cut)),box(x),ref(x))
        else
            return HybridSMC{N,V,T}(SMCg{N,V,T}(cc(x),cv(x),cc_grad(x),cv_grad(x),Intv(x)),box(x),ref(x))
        end
    end
end
