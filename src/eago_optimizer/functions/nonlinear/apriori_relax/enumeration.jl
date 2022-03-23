const ENUM_OUTER_RND = 1E-9
function f_init!(::RelaxMulEnum, g::DAT, b::RelaxCache{V,N,T}) where {V,N,T}
    b.use_apriori_mul = false
    fprop!(Relax(), g, b)
    xp = copy(val(b))
    vlbd = lbd(b)
    vubd = ubd(b)
    xl = copy(vlbd)
    xu = copy(vubd)
    vlbd .-= ENUM_OUTER_RND
    vubd .+= ENUM_OUTER_RND
    for k = 1:node_count(g)
        if !is_num(b, k)
            b._info[k].v = set(b, k)
        end
    end
    s = sparsity(g, 1)
    for i = 0:2^N-1
        s = last(bitstring(i), N)
        for (k,j) in enumerate(sparsity(g, 1))
            b.ic.v.x[j] = s[k] == '1' ? xl[j] : xu[j]
        end
        fprop!(Relax(), g, b)
        for k = node_count(g):-1:1
            if !is_num(b, k)
                b._info[k][i+1] = set(b, k)
            end
        end
    end
    b.ic.v.x .= xp
    vlbd .= xl
    vubd .= xu
    b.use_apriori_mul = true
    fprop!(Relax(), g, b)
    return   
end

_cut_info(t::RelaxMulEnumInner, v, z, x) = z
_cut_info(t::RelaxMulEnum, v, z, x) = z

relax_info(s::RelaxMulEnumInner, n::Int, t::T) where T = MCBoxPnt{2^n,n,T}
relax_info(s::RelaxMulEnum, n::Int, t::T) where T = MCBoxPnt{2^n,n,T}

function estimator_extrema(x::MCBoxPnt{Q,N,T}, y::MCBoxPnt{Q,N,T}, s, dP) where {Q,N,T}
    xmax = maximum(cv, x.box)
    ymax = maximum(cv, y.box) 
    xmin = minimum(cc, x.box)
    ymin = minimum(cc, y.box)
    return xmax, ymax, -xmin, -ymin
end

function estimator_under(xv, yv, x::MCBoxPnt{Q,N,T}, y::MCBoxPnt{Q,N,T}, s, dp, dP, p_rel, p_diam) where {Q,N,T}
    xv.cv, yv.cv, xv.cv_grad, yv.cv_grad
end

function estimator_over(xv, yv, x::MCBoxPnt{Q,N,T}, y::MCBoxPnt{Q,N,T}, s, dp, dP, p_rel, p_diam) where {Q,N,T}
    -xv.cc, -yv.cc, -xv.cc_grad, -yv.cc_grad
end

for d in ALL_ATOM_TYPES
    @eval function fprop!(t::RelaxMulEnum, v::Val{$d}, g::DAT, b::RelaxCache{V,N,T}, k::Int) where {V,N,T<:RelaxTag}
        fprop!(Relax(), v, g, b, k)
    end
end