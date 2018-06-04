Hybrid_Ops1 = [:exp,:exp2,:exp10,:log,:log2,:log10,
               :sqrt,:sqr,:sin,:cos,:inv,:sign,
               :tan,:abs,:acos,:atan,:sinh,
               :cosh,:tanh,:asinh,:atanh,:acosh,:asin]

for k in Hybrid_Ops1
    @eval function ($k)(x::HybridMC)
                   if (hybrid_opts[1].sub_on)
                       return Tighten_Subgrad(HybridMC(($k)(x.SMC)))
                   else
                       return HybridMC(($k)(x.SMC))
                   end
          end
end

+(x::HybridMC) = x
+(x::HybridMC,y::HybridMC) = HybridMC(x.SMC+y.SMC)
+(x::HybridMC,y::T) where {T<:Union{Integer,AbstractFloat}} = HybridMC(x.SMC+y)
+(x::T,y::HybridMC) where {T<:Union{Integer,AbstractFloat}} = HybridMC(x+y.SMC)
-(x::HybridMC) = HybridMC(-x.SMC)

Hybrid_Ops2 = [:-,:*,:/]

for k in Hybrid_Ops2
    @eval ($k)(x::HybridMC,y::T) where {T<:Union{Integer,AbstractFloat}} = HybridMC(($k)(x.SMC,y))
    @eval ($k)(x::T,y::HybridMC) where {T<:Union{Integer,AbstractFloat}} = HybridMC(($k)(x,y.SMC))
    @eval function ($k)(x::HybridMC,y::HybridMC)
             if (hybrid_opts[1].sub_on)
                 return Tighten_Subgrad(HybridMC(($k)(x.SMC,y.SMC)))
             else
                 return HybridMC(($k)(x.SMC,y.SMC))
             end
           end
end

Hybrid_Ops3 = [:min,:max]

for k in Hybrid_Ops3
    @eval function ($k)(x::HybridMC,y::T) where {T<:Union{Integer,AbstractFloat}}
              if (hybrid_opts[1].sub_on)
                  return Tighten_Subgrad(HybridMC(($k)(x.SMC,y)))
              else
                  return HybridMC(($k)(x.SMC,y))
              end
          end
    @eval function ($k)(x::T,y::HybridMC) where {T<:Union{Integer,AbstractFloat}}
              if (hybrid_opts[1].sub_on)
                  return Tighten_Subgrad(HybridMC(($k)(x,y.SMC)))
              else
                   return HybridMC(($k)(x,y.SMC))
              end
          end
    @eval function ($k)(x::HybridMC,y::HybridMC)
              if (hybrid_opts[1].sub_on)
                  return Tighten_Subgrad(HybridMC(($k)(x.SMC,y.SMC)))
              else
                  return HybridMC(($k)(x.SMC,y.SMC))
              end
          end
end

Hybrid_Ops4 = [:pow,:^]
for k in Hybrid_Ops4
    @eval function ($k)(x::HybridMC,y::T) where {T<:Integer}
              if (hybrid_opts[1].sub_on)
                  return Tighten_Subgrad(HybridMC(($k)(x.SMC,y)))
              else
                  return HybridMC(($k)(x.SMC,y))
              end
          end
end

Hybrid_Ops4 = [:pow,:^]
for k in Hybrid_Ops4
    @eval function ($k)(x::HybridMC,y::T) where {T<:AbstractFloat}
              if ~isinteger(y)
                  error("powers must be integral")
              end
              if (hybrid_opts[1].sub_on)
                  return Tighten_Subgrad(HybridMC(($k)(x.SMC,Int(y))))
              else
                  return HybridMC(($k)(x.SMC,Int(y)))
              end
          end
end

one(x::HybridMC) = HybridMC(one(x.SMC))
zero(x::HybridMC) = HybridMC(zero(x.SMC))
dist(x::HybridMC,y::HybridMC) = HybridMC(dist(x.SMC,y.SMC))
real(x::HybridMC) = HybridMC(real(x.SMC))




function convert(::Type{HybridMC{N,V,T}},x::S) where {S<:Integer,N,V,T<:AbstractFloat}
          seed::SVector{N,T} = @SVector zeros(T,N)
          HybridMC{N,V,T}(SMCg{N,V,T}(convert(T,x),convert(T,x),seed,seed,V(convert(V,x)),false))
end
function convert(::Type{HybridMC{N,V,T}},x::S) where {S<:AbstractFloat,N,V,T<:AbstractFloat}
          seed::SVector{N,T} = @SVector zeros(T,N)
          HybridMC{N,V,T}(SMCg{N,V,T}(convert(T,x),convert(T,x),seed,seed,V(convert(V,x)),false))
end
function convert(::Type{HybridMC{N,V,T}},x::S) where {S<:Interval,N,V,T<:AbstractFloat}
          seed::SVector{N,T} = @SVector zeros(T,N)
          HybridMC{N,V,T}(SMCg{N,V,T}(convert(T,x.hi),convert(T,x.lo),seed,seed,convert(V,x),false))
end

promote_rule(::Type{HybridMC{N,V,T}}, ::Type{S}) where {S<:Integer,N,V,T<:AbstractFloat} = HybridMC{N,V,T}
promote_rule(::Type{HybridMC{N,V,T}}, ::Type{S}) where {S<:AbstractFloat,N,V,T<:AbstractFloat} = HybridMC{N,V,T}
promote_rule(::Type{HybridMC{N,V,T}}, ::Type{S}) where {S<:Interval,N,V,T<:AbstractFloat} = HybridMC{N,V,T}
promote_rule(::Type{HybridMC{N,V,T}}, ::Type{S}) where {S<:Real,N,V,T<:AbstractFloat} = HybridMC{N,V,T}
