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

_cut_info(t::RelaxMulEnumInner, v, z, x) = z
_cut_info(t::RelaxMulEnum, v, z, x) = z

relax_info(s::RelaxMulEnumInner, n::Int, t::T) where T = MCBoxPnt{2^n,n,T}
relax_info(s::RelaxMulEnum, n::Int, t::T) where T = MCBoxPnt{2^n,n,T}

const ENUM_OUTER_RND = 1E-9
function f_init!(::RelaxMulEnum, g::DAT, b::RelaxCache{V,N,T}) where {V,N,T}
    b.use_apriori_mul = false
    fprop!(Relax(), g, b)
    xp = copy(b.v.x)
    vlbd = _lbd(b.v)
    vubd = _ubd(b.v)
    xl = copy(vlbd)
    xu = copy(vubd)
    vlbd .-= ENUM_OUTER_RND
    vubd .+= ENUM_OUTER_RND
    for k = 1:_node_count(g)
        if _is_unlocked(b, k) && !_is_num(b, k)
            b._info[k].v = _set(b, k)
        end
    end
    s = _sparsity(g, 1)
    for i = 0:2^N-1
        s = last(bitstring(i), N)
        for (k,j) in enumerate(_sparsity(g, 1))
            b.v.x[j] = s[k] == '1' ? xl[j] : xu[j]
        end
        fprop!(Relax(), g, b)
        for k = _node_count(g):-1:1
            if _is_unlocked(b, k) && !_is_num(b, k)
                b._info[k][i+1] = _set(b, k)
            end
        end
    end
    b.v.x .= xp
    vlbd .= xl
    vubd .= xu
    b.use_apriori_mul = true
    fprop!(RelaxMulEnumInner(), g, b)
    return   
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
        if _is_unlocked(b, k) && !_is_num(b, k)
            b._info[k].v = _set(b, k)
        end
    end
    vlbd .= xl
    vubd .= xu
    fprop!(RelaxMulEnumInner(true), g, b)
    return
end

for F in (EXP, EXP10, LOG, LOG10, POW, MINUS, DIV, PLUS)
    @eval function fprop!(t::RelaxMulEnumInner, v::Val{$F}, g::DAT, b::RelaxCache{V,N,T}, k::Int) where {V,N,T<:RelaxTag}
        fprop!(Relax(), v, g, b, k)
    end
end

function fprop!(t::RelaxMulEnumInner, vt::Variable, g::DAT, b::RelaxCache{V,N,T}, k::Int) where {V,N,T<:RelaxTag}
    i = _first_index(g, k)
    x = _val(b, i)
    l = _lbd(b, i)
    u = _ubd(b, i)
    #@show x, l, u
    z = _var_set(MC{N,T}, _rev_sparsity(g, i, k), x, x, l, u)
    #@show z
    if !_first_eval(b)
        z = z ∩ _interval(b, k)
    end
    #@show z
    _store_set!(b, z, k)
    return
end

function fprop_2!(t::RelaxMulEnumInner, v::Val{MULT}, g::DAT, b::RelaxCache{V,N,T}, k::Int) where {V,N,T<:RelaxTag}
    x = _child(g, 1, k)
    y = _child(g, 2, k)
    x_is_num = _is_num(b, x)
    y_is_num = _is_num(b, y)
    if !x_is_num && y_is_num
        z = _set(b, x)*_num(b, y)
    elseif x_is_num && !y_is_num
        z = _num(b, x)*_set(b, y)
    else
        xs = _set(b, x)
        ys = _set(b, y)
        xinfo = _info(b, x)
        yinfo = _info(b, y)
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
    z = _cut(z, _set(b, k), b.v, b.ϵ_sg, _sparsity(g, k), b.cut, false)
    _store_set!(b, z, k)
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
    for i in _children(g, k)
        if _is_num(b, i)
            znum = znum*_num(b, i)
        else
            xs = _set(b, i)
            xinfo = _info(b, i)
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
    z = _cut(z, _set(b, k), b.v, b.ϵ_sg, _sparsity(g, k), b.cut, false)
    _store_set!(b, z, k)
    return
end

#TODO:
function fprop!(t::RelaxMulEnumInner, v::Val{MULT}, g::DAT, b::RelaxCache{V,N,T}, k::Int) where {V,N,T<:RelaxTag}
    (_arity(g, k) == 2) ? fprop_2!(t, Val(MULT), g, b, k) : fprop_n!(t, Val(MULT), g, b, k)
end