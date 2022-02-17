const ENUM_OUTER_RND = 1E-9
function f_init!(::RelaxMulEnum, g::DAT, b::RelaxCache{V,N,T}) where {V,N,T}
    b.use_apriori_mul = false
    fprop!(Relax(), g, b)
    xp = val(b)
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
    xmin = minimum(cv, x.box)
    ymin = minimum(cv, y.box)
    return xmax, ymax, xmin, ymin
end

function estimator_under(x::MCBoxPnt{Q,N,T}, y::MCBoxPnt{Q,N,T}, s, dp, dP) where {Q,N,T}
    x.v.cv, y.v.cv, x.v.cv_grad, y.v.cv_grad
end

function estimator_over(x::MCBoxPnt{Q,N,T}, y::MCBoxPnt{Q,N,T}, s, dp, dP) where {Q,N,T}
    -x.v.cc, -y.v.cc, -x.v.cc_grad, -y.v.cc_grad
end


#=
function extract_apriori_info(t::RelaxMulEnumInner, x::Vector{MC{N,T}}, Q::Int) where {N,T}
    cv = -Inf
    cc = Inf
    #DEBUG_NL && @show x.box
    for i = 1:Q
        z = x[i]
        cvt = z.cv 
        cct = z.cc
        if cvt > cv
            cv = cvt
        end
        if cct < cc
            cc = cct
        end
    end
    return cv, cc
end
function extract_apriori_info(t::RelaxMulEnumInner, x::MCBoxPnt{Q,N,T}, y::MC{N,T}) where {Q,N,T}
    extract_apriori_info(t, x.box, Q)
end


function fprop!(t::RelaxMulEnum, g::DAT, b::RelaxCache{V,N,T}) where {V,N,T<:RelaxTag}
    vlbd = _lbd(b.v)
    vubd = _ubd(b.v)
    xl = copy(vlbd)
    xu = copy(vubd)
    vlbd .-= ENUM_OUTER_RND
    vubd .+= ENUM_OUTER_RND
    fprop!(RelaxMulEnumInner(false), g, b)
    for k = 1:_node_count(g)
        if !is_num(b, k)
            b._info[k].v = set(b, k)
        end
    end
    vlbd .= xl
    vubd .= xu
    fprop!(RelaxMulEnumInner(true), g, b)
    return
end

fprop!(t::RelaxMulEnumInner, v::Val{EXP}, g::DAT, b::RelaxCache{V,N,T}, k::Int) where {V,N,T<:RelaxTag} = fprop!(Relax(), v, g, b, k)
fprop!(t::RelaxMulEnumInner, v::Val{EXP10}, g::DAT, b::RelaxCache{V,N,T}, k::Int) where {V,N,T<:RelaxTag} = fprop!(Relax(), v, g, b, k)
fprop!(t::RelaxMulEnumInner, v::Val{LOG}, g::DAT, b::RelaxCache{V,N,T}, k::Int) where {V,N,T<:RelaxTag} = fprop!(Relax(), v, g, b, k)
fprop!(t::RelaxMulEnumInner, v::Val{LOG10}, g::DAT, b::RelaxCache{V,N,T}, k::Int) where {V,N,T<:RelaxTag} = fprop!(Relax(), v, g, b, k)
fprop!(t::RelaxMulEnumInner, v::Val{POW}, g::DAT, b::RelaxCache{V,N,T}, k::Int) where {V,N,T<:RelaxTag} = fprop!(Relax(), v, g, b, k)
fprop!(t::RelaxMulEnumInner, v::Val{MINUS}, g::DAT, b::RelaxCache{V,N,T}, k::Int) where {V,N,T<:RelaxTag} = fprop!(Relax(), v, g, b, k)
fprop!(t::RelaxMulEnumInner, v::Val{DIV}, g::DAT, b::RelaxCache{V,N,T}, k::Int) where {V,N,T<:RelaxTag} = fprop!(Relax(), v, g, b, k)
fprop!(t::RelaxMulEnumInner, v::Val{PLUS}, g::DAT, b::RelaxCache{V,N,T}, k::Int) where {V,N,T<:RelaxTag} = fprop!(Relax(), v, g, b, k)
fprop!(t::RelaxMulEnumInner, v::Variable, g::DAT, b::RelaxCache{V,N,T}, k::Int) where {V,N,T<:RelaxTag} = fprop!(Relax(), v, g, b, k) 

function fprop_2!(t::RelaxMulEnumInner, v::Val{MULT}, g::DAT, b::RelaxCache{V,N,T}, k::Int) where {V,N,T<:RelaxTag}
    x = child(g, 1, k)
    y = child(g, 2, k)
    x_is_num = is_num(b, x)
    y_is_num = is_num(b, y)
    if !x_is_num && y_is_num
        z = set(b, x)*num(b, y)
    elseif x_is_num && !y_is_num
        z = num(b, x)*set(b, y)
    else
        xs = set(b, x)
        ys = set(b, y)
        xinfo = info(b, x)
        yinfo = info(b, y)
        xi = t.use_info ? xinfo.v : xs
        yi = t.use_info ? yinfo.v : ys
        xcvU, xccL = extract_apriori_info(t, xinfo, xs)
        ycvU, yccL = extract_apriori_info(t, yinfo, ys)
        xcv = xi.cv; ycv = yi.cv
        xcc = xi.cc; ycc = yi.cc
        xcvg = xi.cv_grad; ycvg = yi.cv_grad
        xccg = xi.cc_grad; yccg = yi.cc_grad
        wIntv = xs.Intv*ys.Intv
        za_l = McCormick.mult_apriori_kernel(xs, ys, wIntv, xcv, ycv, xcvU, ycvU, xcvg, ycvg)
        za_u = McCormick.mult_apriori_kernel(-xs, -ys, wIntv, -xcc, -ycc, -xccL, -yccL, -xccg, -yccg)
        z = (xs*ys) ∩ za_l ∩ za_u
    end
    z = cut(z, set(b, k), b.ic.v, b.ϵ_sg, sparsity(g, k), b.cut, false)
    b[k] = z
    return
end

#TODO: INFO FOR SET BASED ON PASS...
function fprop_n!(t::RelaxMulEnumInner, ::Val{MULT}, g::DAT, b::RelaxCache{V,N,T}, k::Int) where {V,N,T<:RelaxTag}
    ys = zero(MC{N,T})
    znum = one(Float64)
    b._mult_temp.v = one(MC{N,T})
    b._mult_temp.box .= one(MC{N,T})
    count = 0
    first_set = true
    for i in children(g, k)
        if is_num(b, i)
            znum = znum*num(b, i)
        else
            xs = set(b, i)
            xinfo = info(b, i)
            xi = t.use_info ? xinfo.v : xs
            if !first_set
                xcvU, xccL = extract_apriori_info(t, xinfo, xs)
                ycvU, yccL = extract_apriori_info(t, b._mult_temp, ys)
                xcv = xi.cv;       ycv = b._mult_temp.v.cv
                xcc = xi.cc;       ycc = b._mult_temp.v.cc
                xcvg = xi.cv_grad; ycvg = b._mult_temp.v.cv_grad
                xccg = xi.cc_grad; yccg = b._mult_temp.v.cc_grad
                wIntv = xs.Intv*ys.Intv 
                za_l = McCormick.mult_apriori_kernel(xs, ys, wIntv, xcv, ycv, xcvU, ycvU, xcvg, ycvg)
                za_u = McCormick.mult_apriori_kernel(-xs, -ys, wIntv, -xcc, -ycc, -xccL, -yccL, -xccg, -yccg)
                ys = (xs*ys) ∩ za_l ∩ za_u
                b._mult_temp.v *= ys
                b._mult_temp.box .*= xinfo.box
            else 
                first_set = false
                b._mult_temp.v = xi
                b._mult_temp.box .= xinfo.box
                ys = xs
            end
        end
        count += 1
    end
    z = ys*znum
    z = cut(z, set(b, k), b.v, b.ϵ_sg, sparsity(g, k), b.cut, false)
    b[k] = z
    return
end

#TODO:
function fprop!(t::RelaxMulEnumInner, v::Val{MULT}, g::DAT, b::RelaxCache{V,N,T}, k::Int) where {V,N,T<:RelaxTag}
    (arity(g, k) == 2) ? fprop_2!(t, Val(MULT), g, b, k) : fprop_n!(t, Val(MULT), g, b, k)
end
=#