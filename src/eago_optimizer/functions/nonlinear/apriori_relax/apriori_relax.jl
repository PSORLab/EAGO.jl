
include(joinpath(@__DIR__, "affine_arithmetic.jl"))
include(joinpath(@__DIR__, "enumeration.jl"))

function fprop!(t::RelaxMulEnum, vt::Variable, g::DAT, b::RelaxCache{V,N,T}, k::Int) where {V,N,T<:RelaxTag}
    i = _first_index(g, k)
    x = _val(b, i)
    z = _var_set(MC{N,T}, _rev_sparsity(g, i, k), x, x, _lbd(b, i), _ubd(b, i))
    if !_first_eval(b)
        z = z ∩ _interval(b, k)
    end
    _store_set!(b, z, k)    
    if (b.first_eval && b.use_apriori_mul)
        l = [_lbd(b, i) for i = 1:N]
        u = [_ubd(b, i) for i = 1:N]
        for i = 0:2^N-1
            s = bitstring(i)
            v = s[end-n+1:end]
            for j = 1:length(v)
                flag = parse(Bool, v[j]; base = 2)
                if flag
                    MC{N,T}(l[j], Interval(l[j], u[j]), j)
                end
            end
        end
        zs = MCStack{2^N,N,T}(l, u)
        zinfo = MCBoxPnt{2^N,N,T}(z, zs)
        _store_info!(b, zinfo, k)
    end
end

function fprop!(t::RelaxAA, vt::Variable, g::DAT, b::RelaxCache{V,N,T}, k::Int) where {V,N,T<:RelaxTag}
    i = _first_index(g, k)
    x = _val(b, i)
    z = _var_set(MC{N,T}, _rev_sparsity(g, i, k), x, x, _lbd(b, i), _ubd(b, i))
    if !_first_eval(b)
        z = z ∩ _interval(b, k)
    end
    _store_set!(b, z, k)
    (b.first_eval && b.use_apriori_mul) && _store_info!(b, z, k)
end

for RTYPE in (RelaxMulEnum, RelaxAA)

    for (LABEL,f) in ((EXP, :exp), (EXP10, :exp10), (LOG, :log))
        @eval function fprop!(t::$RTYPE, v::Val{$LABEL}, g::DAT, b::RelaxCache{V,N,T}, k::Int) where {V,N,T<:RelaxTag}
            println("blep belp blep = $LABEL")
            if b.first_eval
                xinfo = _info(b, _child(g, 1, k))
                zinfo = ($f)(xinfo)
                _store_info!(b, zinfo, k)
            end
            x = _set(b, _child(g, 1, k))
            z = ($f)(x)
            z = _cut(z, _set(b, k), b.v, zero(Float64), _sparsity(g,k), b.post, b.cut, b.cut_interval)
            zinfo = _info(b, k)
            z = _cut_info(z, zinfo)
            _store_set!(b, z, k)
            return
        end
    end

    @eval function fprop!(t::$RTYPE, v::Val{POW}, g::DAT, b::RelaxCache{V,N,T}, k::Int) where {V,N,T<:RelaxTag}
        println("blep belp blep = POW")
        x = _child(g, 1, k)
        y = _child(g, 2, k)
        x_is_num = _is_num(b, x)
        y_is_num = _is_num(b, y)
        if b.first_eval
            if y_is_num && isone(_num(b, y))
                _store_info!(b, _info(b, x), k)
            elseif y_is_num && iszero(_num(b, y))
                _store_info!(b, one(_info(b, x)), k)
            else
                if !x_is_num && y_is_num
                    zinfo = _info(b, x)^_num(b, y)
                elseif x_is_num && !y_is_num
                    zinfo = _num(b, x)^_info(b, y)
                elseif !x_is_num && !y_is_num
                    zinfo = _info(b, x)^_info(b, y)
                end
                _store_info!(b, zinfo, k)
            end
        end
        if y_is_num && isone(_num(b, y))
            z = _set(b, x)
            _store_set!(b, z, k)
        elseif y_is_num && iszero(_num(b, y))
            _store_set!(b, zero(_set(b, x)), k)
        else
            if !x_is_num && y_is_num
                z = _set(b, x)^_num(b, y)
            elseif x_is_num && !y_is_num
                z = _num(b, x)^_set(b, y)
            elseif !x_is_num && !y_is_num
                z = _set(b, x)^_set(b, y)
            end
            z = _cut(z, _set(b, k), b.v, zero(Float64), _sparsity(g,k), b.post, b.cut, b.cut_interval)
            zinfo = _info(b, k)
            z = _cut_info(z, zinfo)
            _store_set!(b, z, k)
        end
        return
    end

    @eval function fprop!(t::$RTYPE, v::Val{MINUS}, g::DAT, b::RelaxCache{V,N,T}, k::Int) where {V,N,T<:RelaxTag}
        println("blep belp blep = MINUS")
        x = _child(g, 1, k)
        x_is_num = _is_num(b, x)
        if _arity(g, k) == 2
            y = _child(g, 2, k)
            y_is_num = _is_num(b, y)
            if !x_is_num && y_is_num
                z = _set(b, x) - _num(b, y)
            elseif x_is_num && !y_is_num
                z = _num(b, x) - _set(b, y)
            else
                z = _set(b, x) - _set(b, y)
            end
        else
            z = -_set(b, x)
        end
        if b.first_eval
            if _arity(g, k) == 2
                y = _child(g, 2, k)
                y_is_num = _is_num(b, y)
                if !x_is_num && y_is_num
                    zinfo = _info(b, x) - _num(b, y)
                elseif x_is_num && !y_is_num
                    zinfo = _num(b, x) - _info(b, y)
                else
                    zinfo = _info(b, x) - _info(b, y)
                end
            else
                zinfo = -_info(b, x)
            end
            _store_info!(b, zinfo, k)
        end
        zinfo = _info(b, k)
        z = _cut_info(z, zinfo)
        z = _cut(z, _set(b, k), b.v, b.ϵ_sg, _sparsity(g, k), false, b.cut, b.cut_interval)
        _store_set!(b, z, k)
        return
    end

    @eval function fprop!(t::$RTYPE, v::Val{DIV}, g::DAT, b::RelaxCache{V,N,T}, k::Int) where {V,N,T<:RelaxTag}
        println("blep belp blep = DIV")
        x = _child(g, 1, k)
        y = _child(g, 2, k)
        if !_is_num(b, x) && _is_num(b, y)
            z = _set(b, x)/_num(b, y)
        elseif _is_num(b, x) && !_is_num(b, y)
            z = _num(b, x)/_set(b, y)
        else
            z =_set(b, x)/_set(b, y)
        end
        if b.first_eval
            y = _child(g, 2, k)
            y_is_num = _is_num(b, y)
            if !x_is_num && y_is_num
                zinfo = _info(b, x)/_num(b, y)
            elseif x_is_num && !y_is_num
                zinfo = _num(b, x)/_info(b, y)
            else
                zinfo = _info(b, x)/_info(b, y)
            end
            _store_info!(b, zinfo, k)
        end
        zinfo = _info(b, k)
        z = _cut_info(z, zinfo)
        z = _cut(z, _set(b,k), b.v, b.ϵ_sg, _sparsity(g, k), false, b.cut, b.cut_interval)
        _store_set!(b, z, k)
        return
    end

    @eval function fprop_2!(t::$RTYPE, v::Val{PLUS}, g::DAT, b::RelaxCache{V,N,T}, k::Int) where {V,N,T<:RelaxTag}
        println("blep belp blep = PLUS")
        x = _child(g, 1, k)
        y = _child(g, 2, k)
        x_is_num = _is_num(b, x)
        y_is_num = _is_num(b, y)
        if !x_is_num && y_is_num
            z = _set(b, x) + _num(b, y)
        elseif x_is_num && !y_is_num
            z = _num(b, x) + _set(b, y)
        else
            z = _set(b, x) + _set(b, y)
        end
        if b.first_eval
            if !x_is_num && y_is_num
                zinfo = _info(b, x) + _num(b, y)
            elseif x_is_num && !y_is_num
                zinfo = _num(b, x) + _info(b, y)
            else
                zinfo = _info(b, x) + _info(b, y)
            end 
            _store_info!(b, zinfo, k)
        end
        zinfo = _info(b, k)
        z = _cut_info(z, zinfo)
        z = _cut(z, _set(b,k), b.v, b.ϵ_sg, _sparsity(g, k), false, b.cut, b.cut_interval)
        _store_set!(b, z, k)
        return
    end

    @eval function fprop_n!(t::$RTYPE, ::Val{PLUS}, g::DAT, b::RelaxCache{V,N,T}, k::Int) where {V,N,T<:RelaxTag}
        println("blep belp blep = PLUS")
        if b.first_eval
            zinfo = one(MC{N,T})
            znuminfo = one(Float64)
            for i in _children(g, k)
                if _is_num(b, i)
                    znuminfo = znuminfo*_num(b, i)
                else
                    zinfo = zinfo*_info(b, i)
                end
            end
            zinfo = zinfo*znuminfo
            _store_info!(b, zinfo, k)
        end
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
        zinfo = _info(b, k)
        z = _cut_info(z, zinfo)
        z = _cut(z, _set(b, k), b.v, b.ϵ_sg, _sparsity(g, k), false, b.cut, b.cut_interval)
        _store_set!(b, z, k)
        return
    end

    @eval function fprop_2!(t::$RTYPE, v::Val{MULT}, g::DAT, b::RelaxCache{V,N,T}, k::Int) where {V,N,T<:RelaxTag}
        println("blep belp blep = MULT")
        x = _child(g, 1, k)
        y = _child(g, 2, k)
        x_is_num = _is_num(b, x)
        y_is_num = _is_num(b, y)
        if b.first_eval
            if !x_is_num && y_is_num
                zinfo = _info(b, x)*_num(b, y)
            elseif x_is_num && !y_is_num
                zinfo = _num(b, x)*_info(b, y)
            else
                zinfo = _info(b, x)*_info(b, y)
            end 
            _store_info!(b, zinfo, k)
        end
        if !x_is_num && y_is_num
            z = *(_set(b, x), _num(b, y))
        elseif x_is_num && !y_is_num
            z = *(_num(b, x), _set(b, y))
        else
            xv = _set(b, x)
            yv = _set(b, y)
            if b.use_apriori_mul
                xr = _info(b, x)
                yr = _info(b, y)
                xcv, xcvU, xcc, xccL, xcvg, xccg = extract_apriori_info(t, xr, xv)
                ycv, ycvU, ycc, yccL, ycvg, yccg = extract_apriori_info(t, yr, yv) 
                z = mult_apriori_kernel(xv, yv, xv.Intv*yv.Intv, xcv, ycv, xcvU, ycvU, xcc, ycc, xccL, yccL, xcvg, ycvg, xccg, yccg)                                        
            else
                z = xv*yv
            end
        end
        z = _cut(z, _set(b,k), b.v, b.ϵ_sg, _sparsity(g, k), false, b.cut, b.cut_interval)
        _store_set!(b, z, k)
        return
    end

    #TODO:
    function fprop_n!(t::Relax, ::Val{MULT}, g::DAT, b::RelaxCache{V,N,T}, k::Int) where {V,N,T<:RelaxTag}
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
        (b.first_eval && b.use_apriori_mul) && _store_info!(b, z, k)
        return
    end

    for LABEL in (PLUS, MULT)
        @eval function fprop!(t::$RTYPE, v::Val{$LABEL}, g::DAT, b::RelaxCache{V,N,T}, k::Int) where {V,N,T<:RelaxTag}
            println("blep belp blep = $LABEL")
            n = _arity(g, k)
            if n == 2
                return fprop_2!(Relax(), Val($LABEL), g, b, k)
            end
            fprop_n!(Relax(), Val($LABEL), g, b, k)
            return
        end
    end
end