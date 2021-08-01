struct MCStack{Q,N,T}
    v::SVector{Q,MC{N,T}}
end

zero(::Type{MCStack{Q,N,T}}) where {Q,N,T} = MCStack{Q,N,T}(SVector{Q,MC{N,T}}(zeros(MC{N,T},Q)))
zero(x::MCStack{Q,N,T}) where {Q,N,T} = MCStack{Q,N,T}(zero.(x.v))

+(x::MCStack{Q,N,T}, y::MCStack{Q,N,T}) where {Q,N,T} = MCStack{Q,N,T}(x.v + y.v)
+(x::MCStack{Q,N,T}, α::T) where {Q,N,T} = MCStack{Q,N,T}(x.v + α)
+(α::T, x::MCStack{Q,N,T}) where {Q,N,T} = x + α

*(x::MCStack{Q,N,T}, y::MCStack{Q,N,T}) where {Q,N,T} = MCStack{Q,N,T}(x.v .* y.v)
*(x::MCStack{Q,N,T}, α::T) where {Q,N,T} = MCStack{Q,N,T}(x.v * α)
*(α::T, x::MCStack{Q,N,T}) where {Q,N,T} = x*α

-(x::MCStack{Q,N,T}) where {Q,N,T} = MCStack{Q,N,T}(-x.v)
-(x::MCStack{Q,N,T}, α::T) where {Q,N,T} = x + (-α)
-(α::T, x::MCStack{Q,N,T}) where {Q,N,T} = α + (-x)

/(x::MCStack{Q,N,T}, α::T) where {Q,N,T} = MCStack{Q,N,T}(x.v/α)
/(α::T, x::MCStack{Q,N,T}) where {Q,N,T} = MCStack{Q,N,T}(α*inv.(x.v))

^(x::MCStack{Q,N,T}, n::Integer) where {Q,N,T} = MCStack{Q,N,T}(x.v.^n)
inv(x::MCStack{Q,N,T}) where {Q,N,T} = MCStack{Q,N,T}(inv.(x.v))

log(x::MCStack{Q,N,T}) where {Q,N,T} = MCStack{Q,N,T}(log.(x.v))
log10(x::MCStack{Q,N,T}) where {Q,N,T} = MCStack{Q,N,T}(log10.(x.v))

exp(x::MCStack{Q,N,T}) where {Q,N,T} = MCStack{Q,N,T}(exp.(x.v))
exp10(x::MCStack{Q,N,T}) where {Q,N,T} = MCStack{Q,N,T}(exp10.(x.v))

struct MCBoxPnt{Q,N,T}
    v::MC{N,T}
    box::MCStack{Q,N,T}
end

zero(::Type{MCBoxPnt{Q,N,T}}) where {Q,N,T} = MCBoxPnt{Q,N,T}(zero(MC{N,T}), zero(MCStack{Q,N,T}))
zero(x::MCBoxPnt{Q,N,T}) where {Q,N,T} = MCBoxPnt{Q,N,T}(zero(x.v), zero(x.box))

+(x::MCBoxPnt{Q,N,T}, y::MCBoxPnt{Q,N,T}) where {Q,N,T} = MCBoxPnt{Q,N,T}(x.v + y.v, x.box + y.box)
+(x::MCBoxPnt{Q,N,T}, α::T) where {Q,N,T} = MCBoxPnt{Q,N,T}(x.v + α, x.box + α)
+(α::T, x::MCBoxPnt{Q,N,T}) where {Q,N,T} = x + α

*(x::MCBoxPnt{Q,N,T}, y::MCBoxPnt{Q,N,T}) where {Q,N,T} = MCBoxPnt{Q,N,T}(x.v*y.v, x.box*y.box)
*(x::MCBoxPnt{Q,N,T}, α::T) where {Q,N,T} = MCBoxPnt{Q,N,T}(x.v * α, x.box * α)
*(α::T, x::MCBoxPnt{Q,N,T}) where {Q,N,T} = x*α

-(x::MCBoxPnt{Q,N,T}) where {Q,N,T} = MCBoxPnt{Q,N,T}(-x.v, -x.box)
-(x::MCBoxPnt{Q,N,T}, α::T) where {Q,N,T} = x + (-α)
-(α::T, x::MCBoxPnt{Q,N,T}) where {Q,N,T} = α + (-x)

/(x::MCBoxPnt{Q,N,T}, α::T) where {Q,N,T} = MCBoxPnt{Q,N,T}(x.v/α, x.box/α)
/(α::T, x::MCBoxPnt{Q,N,T}) where {Q,N,T} = MCBoxPnt{Q,N,T}(α*inv(x.v), α*inv(x.box))

^(x::MCBoxPnt{Q,N,T}, n::Integer) where {Q,N,T} = MCBoxPnt{Q,N,T}(x.v^n, x.box^n)
inv(x::MCBoxPnt{Q,N,T}) where {Q,N,T} = MCBoxPnt{Q,N,T}(inv(x.v), inv(x.box))

log(x::MCBoxPnt{Q,N,T}) where {Q,N,T} = MCBoxPnt{Q,N,T}(log(x.v), log(x.box))
log10(x::MCBoxPnt{Q,N,T}) where {Q,N,T} = MCBoxPnt{Q,N,T}(log10(x.v), log10(x.box))

exp(x::MCBoxPnt{Q,N,T}) where {Q,N,T} = MCBoxPnt{Q,N,T}(exp(x.v), exp(x.box))
exp10(x::MCBoxPnt{Q,N,T}) where {Q,N,T} = MCBoxPnt{Q,N,T}(exp10(x.v), exp10(x.box))

function extremal_mc(x::MCStack{Q,N,T}) where {Q,N,T}
    cv = -Inf
    cc = Inf
    min_i = -1
    max_i = -1
    for i = 1:Q
        z = x.v[i]
        cvt = z.cv 
        cct = z.cc
        if cvt > cv
            min_i = i
            cv = cvt
        end
        if cct < cc
            max_i = i
            cc = cct
        end
    end
    cv_grad = x.v[min_i].cv_grad
    cc_grad = x.v[max_i].cc_grad
    Intv = x.v[1].Intv
    cnst = x.v[1].cnst
    return MC{N,T}(cv, cc, Intv, cv_grad, cc_grad, cnst)
end

function extract_apriori_info(t::RelaxMulEnum, x::MCStack{Q,N,T}, y::MC{N,T}) where {Q,N,T}
    z = extremal_mc(x::MCStack{Q,N,T})
    return y.cv, z.cv, y.cc, z.cc, y.cv_grad, y.cc_grad
end

relax_info(s::RelaxMulEnum, n::Int, t::T) where T = MCBoxPnt{2^n,n,T}
f_init!(::RelaxMulEnum, g::DAT, b::RelaxCache) = nothing