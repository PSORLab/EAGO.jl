# Defines type
"""
"""
struct HybridMC{N,V,T}
    SMC::SMCg{N,V,T}
    box::SVector{V}
    ref::SVector{T}
end

# Access functions
"""
"""
cc(x::HybridMC) = x.SMC.cc

"""
"""
cv(x::HybridMC) = x.SMC.cv

"""
"""
cc_grad(x::HybridMC) = x.SMC.cc_grad

"""
"""
cv_grad(x::HybridMC) = x.SMC.cv_grad

"""
"""
Intv(x::HybridMC) = x.SMC.Intv

"""
"""
box(x::HybridMC) = x.box

"""
"""
ref(x::HybridMC) = x.reg

boxlo(x.HybridMC) = lo(x.box)
boxhi(x.HybridMC) = hi(x.box)

# constructors

# conversion/promotion functions
