
function affine_expand_del(dx::Vector{Float64}, fx0::Float64, ∇fx0::SVector{N,Float64}, s::Vector{Int}) where N
    v = fx0
    for i=1:N
        v += ∇fx0[i]*dx[s[i]]
    end
    return v
end
function affine_expand_del(dx::Vector{Interval{Float64}}, fx0::Float64, ∇fx0::SVector{N,Float64}, s::Vector{Int}) where N
    v = fx0
    for i = 1:N
        t = ∇fx0[i]
        tdx = dx[s[i]]
        v += t > 0.0 ? t*tdx.hi : t*tdx.lo
    end
    return v
end

function affine_expand(x::Vector{Float64}, x0::Vector{Float64}, fx0::Float64, ∇fx0::SVector{N,Float64}) where N
    v = fx0
    for i=1:N
        v += ∇fx0[i]*(x[i] - x0[i])
    end
    return v
end
function affine_expand(x::Vector{Interval{Float64}}, x0::Vector{Float64}, fx0::Float64, ∇fx0::SVector{N,Float64}) where N
    v = Interval{Float64}(fx0)
    for i=1:N
        v += ∇fx0[i]*(x[i] - x0[i])
    end
    return v
end

function expand_set(::Type{MC{N2,T}}, x::MC{N1,T}, fsparse::Vector{Int}, subsparse::Vector{Int}, cv_buffer::Vector{Float64}, cc_buffer::Vector{Float64}) where {N1, N2, T<:RelaxTag}
    cvg = x.cv_grad
    ccg = x.cc_grad
    xcount = 1
    xcurrent = subsparse[1]
    for i = 1:N2
        if fsparse[i] === xcurrent
            cv_buffer[i] = cvg[xcount]
            cc_buffer[i] = ccg[xcount]
            xcount += 1
            if xcount <= N1
            xcurrent = subsparse[xcount]
            else
                break
            end
        else
            cv_buffer[i] = zero(Float64)
        end
    end
    cv_grad = SVector{N2,Float64}(cv_buffer)
    cc_grad = SVector{N2,Float64}(cc_buffer)
    return MC{N2,T}(x.cv, x.cc, x.Intv, cv_grad, cc_grad, x.cnst)
end

function set_value_post(z::MC{N,T}, v::VariableValues{Float64}, s::Vector{Int}, ϵ::Float64) where {V,N,T<:RelaxTag}
    l = z.cv
    u = z.cc
    lower_refinement = true
    upper_refinement = true
    @inbounds for i = 1:N
        cv_val = z.cv_grad[i]
        cc_val = z.cc_grad[i]
        i_sol = s[i]
        x_z = v.x[i_sol]
        lower_bound = _lbd(v, i_sol)
        upper_bound = _ubd(v, i_sol)
        if lower_refinement
            if cv_val > zero(Float64)
                if isinf(lower_bound)
                    !upper_refinement && break
                    lower_refinement = false
                else
                    l += cv_val*(lower_bound - x_z)
                end
            else
                if isinf(upper_bound)
                    !upper_refinement && break
                    lower_refinement = false
                else
                    l += cv_val*(upper_bound - x_z)
                end
            end
        end
        if upper_refinement
            if cc_val > zero(Float64)
                if isinf(lower_bound)
                    !lower_refinement && break
                    upper_refinement = false
                else
                    u += cc_val*(upper_bound - x_z)
                end
            else
                if isinf(upper_bound)
                    !lower_refinement && break
                    upper_refinement = false
                else
                    u += cc_val*(lower_bound - x_z)
                end
            end
        end
    end

    if lower_refinement && (z.Intv.lo + ϵ > l)
        l = z.Intv.lo
    elseif !lower_refinement
        l = z.Intv.lo
    else
        l -= ϵ
    end

    if upper_refinement && (z.Intv.hi - ϵ < u)
        u = z.Intv.hi
    elseif !upper_refinement
        u = z.Intv.hi
    else
        u += ϵ
    end

    return MC{N,T}(z.cv, z.cc, Interval{Float64}(l, u), z.cv_grad, z.cc_grad, z.cnst)
end

"""
$(FUNCTIONNAME)

Intersects the new set valued operator with the prior and performs affine bound tightening

- First forward pass: `post` should be set by user option, `is_intersect` should be false
  so that the tape overwrites existing values, and the `interval_intersect` flag could be set
  to either value.
- Forward CP pass (assumes same reference point): `post` should be set by user option,
  `is_intersect` should be true so that the tape intersects with  existing values, and the
  `interval_intersect` flag should be false.
- Forward CP pass (assumes same reference point): `post` should be set by user option,
  `is_intersect` should be true so that the tape intersects with existing values, and the
  `interval_intersect` flag should be false.
- Subsequent forward passes at new points: post` should be set by user option,
  `is_intersect` should be true so that the tape intersects with existing values, and the
  `interval_intersect` flag should be `true` as predetermined interval bounds are valid but
   the prior values may correspond to different points of evaluation.
"""
function cut(x::MC{N,T}, z::MC{N,T}, v::VariableValues, ϵ::Float64, s::Vector{Int}, cflag::Bool, pflag::Bool) where {N,T<:RelaxTag}
    (pflag & cflag)  && (return set_value_post(x ∩ z.Intv, v, s, ϵ))
    (pflag & !cflag) && (return set_value_post(x, v, s, ϵ))
    (pflag & cflag)  && (return x ∩ z.Intv)
    return x
end