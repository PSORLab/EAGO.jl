"""
$(FUNCTIONNAME)

An operator that cuts the `x_mc` object using the `x_mc_int bounds` in a
differentiable or nonsmooth fashion to achieve a composite relaxation within
`x_mc_int`.
"""
function final_cut(x_mc::MC{N,NS}, x_mc_int::MC{N,NS}) where {N}
    Intv = x_mc.Intv ∩ x_mc_int.Intv
    if (x_mc.cc < x_mc_int.cc)
      cc = x_mc.cc
      cc_grad::SVector{N,Float64} = x_mc.cc_grad
    else
      cc = x_mc_int.cc
      cc_grad = x_mc_int.cc_grad
    end
    if (x_mc.cv > x_mc_int.cv)
      cv = x_mc.cv
      cv_grad::SVector{N,Float64} = x_mc.cv_grad
    else
      cv = x_mc_int.cv
      cv_grad = x_mc_int.cv_grad
    end
    x_out::MC{N,NS} = MC{N,NS}(cv, cc,(x_mc.Intv ∩ x_mc_int.Intv), cv_grad, cc_grad, x_mc.cnst)
  return x_out
end
function final_cut(x_mc::MC{N,Diff},x_mc_int::MC{N,Diff}) where {N}
  x_out = smooth_cut(x_mc,x_mc_int)
  return x_out
end

abstract type AbstractContractorMC end
struct NewtonGS <: AbstractContractorMC end
struct KrawczykCW <: AbstractContractorMC end

abstract type AbstractPreconditionerMC end
function preconditioner_storage(x::AbstractPreconditionerMC, t::T) where T<:RelaxTag
    error("Must define function that generates appropriate storage type for preconditioner")
end

abstract type AbstractMCCallback end
struct MCCallback{FH <: Function, FJ <: Function, C <: AbstractContractorMC,
                  PRE <: AbstractPreconditionerMC, N, T<:RelaxTag,
                  AMAT <: AbstractMatrix} <:  AbstractMCCallback
    h!::FH
    hj!::FJ
    H::Vector{MC{N,T}}
    J::AMAT
    xmid::Vector{Float64}
    X::Vector{Interval{Float64}}
    P::Vector{Interval{Float64}}
    nx::Int
    np::Int
    lambda::Float64
    eps::Float64
    kmax::Int
    pref_mc::Vector{MC{N,T}}
    p_mc::Vector{MC{N,T}}
    x0_mc::Vector{MC{N,T}}
    x_mc::Vector{MC{N,T}}
    xa_mc::Vector{MC{N,T}}
    xA_mc::Vector{MC{N,T}}
    aff_mc::Vector{MC{N,T}}
    z_mc::Vector{MC{N,T}}
    contractor::C
    preconditioner::PRE
    apply_precond::Bool
    param::Vector{Vector{MC{N,T}}}
end
function MCCallback(h!::FH, hj!::FJ, nx::Int, np::Int,
                    contractor::S = NewtonGS(),
                    preconditioner::T = DenseMidInv(zeros(Float64,1,1), zeros(Interval{Float64},1), 1, 1),
                    relax_tag::TAG = NS()) where {FH <: Function,
                                                  FJ <: Function,
                                                  S <: AbstractContractorMC,
                                                  T,
                                                  TAG <: RelaxTag}
    H = zeros(MC{np,TAG}, (nx,))
    xmid = zeros(Float64, (nx,))
    P = zeros(Interval{Float64}, (np,))
    X = zeros(Interval{Float64}, (nx,))
    lambda = 0.5
    eps = 0.0
    kmax = 2
    p_ref = zeros(MC{np,TAG}, (np,))
    p_mc = zeros(MC{np,TAG}, (np,))
    x0_mc = zeros(MC{np,TAG}, (nx,))
    x_mc = zeros(MC{np,TAG}, (nx,))
    xa_mc = zeros(MC{np,TAG}, (nx,))
    xA_mc = zeros(MC{np,TAG}, (nx,))
    aff_mc = zeros(MC{np,TAG}, (nx,))
    z_mc = zeros(MC{np,TAG}, (nx,))
    contractor = NewtonGS()
    preconditioner = preconditioner(h!, hj!, nx, np)
    J = preconditioner_storage(preconditioner, relax_tag)
    apply_precond = true
    param = fill(zeros(MC{np,TAG}, (nx, )), (kmax,))

    return MCCallback{FH, FJ, NewtonGS, DenseMidInv, np, TAG, typeof(J)}(h!, hj!, H, J, xmid, X, P, nx, np,
                                                                         lambda, eps, kmax, p_ref, p_mc, x0_mc, x_mc,
                                                                         xa_mc, xA_mc, aff_mc, z_mc,
                                                                         contractor, preconditioner,
                                                                         apply_precond, param)
end
function (d::MCCallback)()
    d.h!(d.H, d.z_mc, d.p_mc)
    d.hj!(d.J, d.aff_mc, d.p_mc)
    return
end
# performance optimized

include("preconditioner/dense.jl")

"""
$(FUNCTIONNAME)

Computates the affine relaxations of the state variable.
"""
function affine_exp!(x::S, p::Vector{MC{N,T}}, d::MCCallback) where {S, N, T<:RelaxTag}

    S1 = zero(MC{N,T})
    S2 = zero(MC{N,T})
    S3 = zero(MC{N,T})
    @inbounds for i = 1:d.nx
        S1 = zero(MC{N,T})
        S2 = zero(MC{N,T})
        S3 = zero(MC{N,T})
        @inbounds for j = 1:N
            S1 += (p[j]-d.pref_mc[j])*x[i].cv_grad[j]
            S2 += (p[j]-d.pref_mc[j])*x[i].cc_grad[j]
            S3 += (d.lambda*x[i].cv_grad[j]+(1.0-d.lambda)*x[i].cc_grad[j])*(p[j]-d.pref_mc[j])
        end
        temp1 = x[i].cv + S1
        temp2 = x[i].cc + S2
        temp3 = d.lambda*x[i].cv+(1.0-d.lambda)*x[i].cc+S3
        d.xa_mc[i] = MC{N,T}(temp1.cv, temp1.cc, temp1.Intv, x[i].cv_grad, x[i].cv_grad, S1.cnst)
        d.xA_mc[i] = MC{N,T}(temp2.cv, temp2.cc, temp2.Intv, x[i].cc_grad, x[i].cc_grad, S2.cnst)
        d.z_mc[i] = MC{N,T}(temp3.cv, temp3.cc, temp3.Intv,
                            d.lambda*x[i].cv_grad+(1.0-d.lambda)*x[i].cc_grad,
                            d.lambda*x[i].cv_grad+(1.0-d.lambda)*x[i].cc_grad, S3.cnst)
    end
    return
end
# affine_exp! currently optimized....

"""
$(FUNCTIONNAME)

Corrects the relaxation of the state variable `x_mc` if the affine relaxation,
"""
function correct_exp!(d::MCCallback{FH,FJ,C,PRE,N,T}) where {FH <: Function,
                                                             FJ <: Function,
                                                             C, PRE, N,
                                                             T<:RelaxTag}
    zero_grad = zeros(SVector{N,Float64})
    @inbounds for i = 1:d.nx
        if (d.z_mc[i].Intv.lo - d.eps < d.X[i].lo) && (d.z_mc[i].Intv.hi + d.eps > d.X[i].hi)
            d.x_mc[i] = MC{N,T}(d.X[i])
        elseif (d.z_mc[i].Intv.lo - d.eps < d.X[i].lo)
            d.x_mc[i] = MC{N,T}(d.X[i].lo, d.x_mc[i].cc, Interval(d.X[i].lo, d.x_mc[i].Intv.hi),
                                zero_grad, d.x_mc[i].cc_grad, d.x_mc[i].cnst)
        elseif (d.z_mc[i].Intv.hi + d.eps > d.X[i].hi)
            d.x_mc[i] = MC{N,T}(d.x_mc[i].cv, d.X[i].hi, Interval(d.x_mc[i].Intv.lo, d.X[i].hi),
                                d.x_mc[i].cv_grad, zero_grad, d.x_mc[i].cnst)
        end
    end
    return
end
# correct_exp! currently optimized....

include("contract.jl")

"""
$(FUNCTIONNAME)
"""
function precond_and_contract!(callback!::MCCallback{FH,FJ,C,PRE,N,T}) where {FH <: Function,
                                                                              FJ <: Function,
                                                                              C <: AbstractContractorMC,
                                                                              PRE <: AbstractPreconditionerMC,
                                                                              N, T<:RelaxTag}
    @. callback!.aff_mc = MC{N,T}(cv(callback!.xa_mc), cc(callback!.xA_mc))
    callback!()
    if callback!.apply_precond
        precondition!(callback!.preconditioner, callback!.H, callback!.J)
    end
    contract!(callback!.contractor, callback!)
    return
end
# precond_and_contract! currently optimized....

"""
$(FUNCTIONNAME)

Populates `x_mc`, `xa_mc`, `xA_mc`, and `z_mc` with affine bounds.
"""
function populate_affine!(d::MCCallback{FH,FJ,C,PRE,N,T}, interval_bnds::Bool) where {FH <: Function,
                                                                                      FJ <: Function,
                                                                                      C <: AbstractContractorMC,
                                                                                      PRE <: AbstractPreconditionerMC,
                                                                                      N, T<:RelaxTag}
    if interval_bnds
        @inbounds for i in 1:d.nx
            d.x_mc[i] = MC{N,T}(d.X[i].lo, d.X[i].hi)
            d.xa_mc[i] = MC{N,T}(d.X[i].lo, d.X[i].lo)
            d.xA_mc[i] = MC{N,T}(d.X[i].hi, d.X[i].hi)
            d.z_mc[i] = d.lambda*d.xa_mc[i] + (1.0 - d.lambda)*d.xA_mc[i]
         end
    end
    return
end
# populate_affine! currently optimized....

"""
$(FUNCTIONNAME)
"""
function gen_expansion_params!(d::MCCallback, interval_bnds::Bool = true) where {N, T<:RelaxTag}
    populate_affine!(d, interval_bnds)
    @. d.param[1] = d.x_mc
    for k = 2:d.kmax
        precond_and_contract!(d)
        affine_exp!(d.x_mc, d.pref_mc, d)
        correct_exp!(d)
        @. d.param[k] = d.x_mc
    end
    return
end

"""
$(FUNCTIONNAME)
"""
function implicit_relax_h!(d::MCCallback, interval_bnds::Bool = true) where {N, T<:RelaxTag}
    populate_affine!(d, interval_bnds)
    for k=1:d.kmax
        affine_exp!(d.param[k], d.p_mc, d)
        precond_and_contract!(d)
    end
    return
end
