mutable struct MCBoxPnt{Q,N,T}
    v::MC{N,T}
    box::Vector{MC{N,T}}
end
function zero(::Type{MCBoxPnt{Q,N,T}}) where {Q,N,T}
    MCBoxPnt{Q,N,T}(zero(MC{N,T}), zeros(MC{N,T}, Q))
end
function setindex!(d::MCBoxPnt{Q,N,T}, x::MC{N,T}, i::Int) where {Q,N,T}
    d.box[i] = x
end

function extremal_mc(x::MCBoxPnt{Q,N,T}) where {Q,N,T}
    cv = -Inf
    cc = Inf
    min_i = -1
    max_i = -1
    for i = 1:Q
        z = x.box[i]
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
    cv_grad = x.box[min_i].cv_grad
    cc_grad = x.box[max_i].cc_grad
    Intv = x.v.Intv
    cnst = x.v.cnst
    return MC{N,T}(cv, cc, Intv, cv_grad, cc_grad, cnst)
end

function extract_apriori_info(t::RelaxMulEnum, x::MCBoxPnt{Q,N,T}, y::MC{N,T}) where {Q,N,T}
    z = extremal_mc(x)
    return y.cv, z.cv, y.cc, z.cc, y.cv_grad, y.cc_grad
end
_cut_info(t::RelaxMulEnum, v, z, x) = z

relax_info(s::RelaxMulEnum, n::Int, t::T) where T = MCBoxPnt{2^n,n,T}

function f_init!(::RelaxMulEnum, g::DAT, b::RelaxCache{V,N,T}) where {V,N,T}
    fprop!(Relax(), g, b)
    for k = 1:_node_count(g)
        if _is_unlocked(b, k) && !_is_num(b, k)
            b._info[k].v = _set(b, k)
        end
    end
    xp = copy(b.v.x) 
    s = _sparsity(g, 1)
    for i = 0:2^N-1
        s = last(bitstring(i), N)
        for j in _sparsity(g, 1)
            b.v.x[j] = s[j] == '1' ? _lbd(b, j) : _ubd(b, j)
        end
        fprop!(Relax(), g, b)
        for k = 1:_node_count(g)
            if _is_unlocked(b, k) && !_is_num(b, k)
                b._info[k][i+1] = _set(b, k)
            end
        end
    end
    b.v.x .= xp
    return   
end

for F in (EXP, EXP10, LOG, LOG10, POW, MINUS, DIV, PLUS)
    @eval fprop!(t::RelaxMulEnum, v::Val{$F}, g::DAT, b::RelaxCache{V,N,T}, k::Int) where {V,N,T<:RelaxTag} = fprop!(Relax(), v, g, b, k)
end

function fprop!(t::RelaxMulEnum, vt::Variable, g::DAT, b::RelaxCache{V,N,T}, k::Int) where {V,N,T<:RelaxTag}
    i = _first_index(g, k)
    x = _val(b, i)
    z = _var_set(MC{N,T}, _rev_sparsity(g, i, k), x, x, _lbd(b, i), _ubd(b, i))
    if !_first_eval(b)
        z = z ∩ _interval(b, k)
    end
    _store_set!(b, z, k)
    return
end

function fprop_2!(t::RelaxMulEnum, v::Val{MULT}, g::DAT, b::RelaxCache{V,N,T}, k::Int) where {V,N,T<:RelaxTag}
    x = _child(g, 1, k)
    y = _child(g, 2, k)
    x_is_num = _is_num(b, x)
    y_is_num = _is_num(b, y)
    if !x_is_num && y_is_num
        z = *(_set(b, x), _num(b, y))
    elseif x_is_num && !y_is_num
        z = *(_num(b, x), _set(b, y))
    else
        xv = _set(b, x)
        yv = _set(b, y)
        xcv, xcvU, xcc, xccL, xcvg, xccg = extract_apriori_info(t, _info(b, x), xv)
        ycv, ycvU, ycc, yccL, ycvg, yccg = extract_apriori_info(t, _info(b, y), yv) 
        z = mult_apriori_kernel(xv, yv, xv.Intv*yv.Intv, xcv, ycv, xcvU, ycvU, xcc, ycc, xccL, yccL, xcvg, ycvg, xccg, yccg)                                        
    end
    z = _cut(z, _set(b,k), b.v, b.ϵ_sg, _sparsity(g, k), false, b.cut, b.cut_interval)
    _store_set!(b, z, k)
    return
end

#TODO:
function fprop_n!(t::RelaxMulEnum, ::Val{MULT}, g::DAT, b::RelaxCache{V,N,T}, k::Int) where {V,N,T<:RelaxTag}
    println("blep belp blep = MULT")
    z = one(MC{N,T})
    znum = one(Float64)
    count = 0
    for i in _children(g, k)
        if _is_num(b, i)
            znum = znum*_num(b, i)
        else
            z = z*_set(b, i)
        end
        count += 1
    end
    z = z*znum
    z = _cut(z, _set(b, k), b.v, b.ϵ_sg, _sparsity(g, k), false, b.cut, b.cut_interval)
    _store_set!(b, z, k)
    return
end

function fprop!(t::RelaxMulEnum, v::Val{MULT}, g::DAT, b::RelaxCache{V,N,T}, k::Int) where {V,N,T<:RelaxTag}
    (_arity(g, k) == 2) ? fprop_2!(t, Val(MULT), g, b, k) : fprop_n!(t, Val(MULT), g, b, k)
end