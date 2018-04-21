# Defines type
"""

Defines the hybridMC type used for constructing nonstandard McCormick relaxations
"""
struct HybridMC{N,V,T}
    SMC::SMCg{N,V,T}
end

############################## constructors ####################################

#################### conversion/promotion functions ############################

###########################  Access functions  #################################
"""
    cc
Access concave relaxation of hybrid McCormick object.
"""
cc(x::HybridMC) = x.SMC.cc

"""
    cv

Access convex relaxation of hybrid McCormick object.
"""
cv(x::HybridMC) = x.SMC.cv

"""
    cc_grad
Access concave (sub)gradient of hybrid McCormick object.
"""
cc_grad(x::HybridMC) = x.SMC.cc_grad

"""
    cv_grad
Access convex (sub)gradient of hybrid McCormick object.
"""
cv_grad(x::HybridMC) = x.SMC.cv_grad

"""
    Intv
Access the interval bounds of the hybrid McCormick object.
"""
Intv(x::HybridMC) = x.SMC.Intv
