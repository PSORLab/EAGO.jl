# Copyright (c) 2018: Matthew Wilhelm & Matthew Stuber.
# This code is licensed under MIT license (see LICENSE.md for full details)
#############################################################################
# EAGO
# A development environment for robust and global optimization
# See https://github.com/PSORLab/EAGO.jl
#############################################################################
# src/eago_optimizer/evaluator/passes.jl
# Functions used to compute forward pass of nonlinear functions which include:
# set_value_post, overwrite_or_intersect, forward_pass_kernel, associated blocks
#############################################################################

function _var_set(::Type{MC{N,T}}, i::Int, x_cv::Float64, x_cc::Float64, l::Float64, u::Float64) where {V,N,T<:RelaxTag}
    v = seed_gradient(i, Val(N))
    return MC{N,T}(x_cv, x_cc, Interval{Float64}(l, u), v, v, false)
end
function fprop!(t::Relax, vt::Variable, g::DAT, b::RelaxCache{V,N,T}, k::Int) where {V,N,T<:RelaxTag}
    i = _first_index(g, k)
    x = _val(b, i)
    z = _var_set(MC{N,T}, _rev_sparsity(g, i, k), x, x, _lbd(b, i), _ubd(b, i))
    if !_first_eval(b)
        z = z ∩ _interval(b, k)
    end
    _store_set!(b, z, k)
end

function fprop!(t::Relax, ex::Subexpression, g::DAT, c::RelaxCache{V,N,T}, k::Int) where {V,N,T<:RelaxTag}
    z = _subexpression_set(c, _first_index(g, k))
    _store_set!(c, z, k)
end

for (f, F, fc) in ((:fprop_2!, PLUS, :+), (:fprop_2!, MIN, :min), (:fprop_2!, MAX, :max), (:fprop!, DIV, :/))
    eval(quote
        function ($f)(t::Relax, v::Val{$F}, g::DAT, b::RelaxCache{V,N,T}, k::Int) where {V,N,T<:RelaxTag}
            x = _child(g, 1, k)
            y = _child(g, 2, k)
            if !_is_num(b, x) && _is_num(b, y)
                z = ($fc)(_set(b, x), _num(b, y))
            elseif _is_num(b, x) && !_is_num(b, y)
                z = ($fc)(_num(b, x), _set(b, y))
            else
                z = ($fc)(_set(b, x), _set(b, y))
            end
            z = _cut(z, _set(b,k), b.v, b.ϵ_sg, _sparsity(g, k), b.cut, false)
            _store_set!(b, z, k)
        end
    end)
end
function fprop!(t::Relax, v::Val{MINUS}, g::DAT, b::RelaxCache{V,N,T}, k::Int) where {V,N,T<:RelaxTag}
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
    z = _cut(z, _set(b, k), b.v, b.ϵ_sg, _sparsity(g, k), b.cut, false)
    _store_set!(b, z, k)
end


function fprop_2!(t::Relax, v::Val{MULT}, g::DAT, b::RelaxCache{V,N,T}, k::Int) where {V,N,T<:RelaxTag}

    x = _child(g, 1, k)
    y = _child(g, 2, k)
    node_x = g.nodes[x]
    node_y = g.nodes[y]

    x_is_num = _is_num(b, x)
    y_is_num = _is_num(b, y)
    if !x_is_num && y_is_num
        z = *(_set(b, x), _num(b, y))
    elseif x_is_num && !y_is_num
        z = *(_num(b, x), _set(b, y))
    else
        xv = _set(b, x)
        yv = _set(b, y)
        if b.use_apriori_mul
            dp = b.dp
            dP = b.dP
            s = _sparsity(g, 1)

            xr = _info(b, x)
            yr = _info(b, y)

            xr_cv = cv(xr)                   # GOOD
            xr_cc = cc(xr)                   # GOOD
            yr_cv = cv(yr)                   # GOOD
            yr_cc = cc(yr)                   # GOOD
            xr_cv_grad = cv_grad(xr)         # GOOD
            xr_cc_grad = cc_grad(xr)         # GOOD
            yr_cv_grad = cv_grad(yr)         # GOOD
            yr_cc_grad = cc_grad(yr)         # GOOD

            xrn_cc = -xr_cc
            yrn_cc = -yr_cc
            xrn_cc_grad = -xr_cc_grad
            yrn_cc_grad = -yr_cc_grad

            t3 = affine_expand_del(dP, xr_cv, xr_cv_grad, s)
            t4 = affine_expand_del(dP, yr_cv, yr_cv_grad, s)
            s3 = affine_expand_del(dP, xrn_cc, xrn_cc_grad, s)
            s4 = affine_expand_del(dP, yrn_cc, yrn_cc_grad, s)

            z = xv*yv
            wIntv = z.Intv
            if (t3 < x.Intv.hi) || (t4 < z.Intv.hi)
                t1 = affine_expand_del(dp, xr_cv, xr_cv_grad, s)
                t2 = affine_expand_del(dp, yr_cv, yr_cv_grad, s)
                za_l = McCormick.mult_apriori_kernel(xv, yv, wIntv, t1, t2, t3, t4, xr_cv_grad, yr_cv_grad)
                z = z ∩ za_l
            end
            if (s3 < -x.Intv.lo) || (s4 < -z.Intv.lo)
                s1 = affine_expand_del(dp, xrn_cc, xrn_cc_grad, s)
                s2 = affine_expand_del(dp, yrn_cc, yrn_cc_grad, s)
                za_u = McCormick.mult_apriori_kernel(-xv, -yv, wIntv, s1, s2, s3, s4, xrn_cc_grad, yrn_cc_grad)
                z = z ∩ za_u
            end
        else
            z = xv*yv
        end
    end
    z = _cut(z, _set(b,k), b.v, b.ϵ_sg, _sparsity(g, k), b.cut, false)
    _store_set!(b, z, k)
end

function fprop_n!(t::Relax, v::Val{PLUS}, g::DAT, b::RelaxCache{V,N,T}, k::Int) where {V,N,T<:RelaxTag}
    z = zero(MC{N,T})
    znum = 0.0
    for i in _children(g, k)
        if _is_num(b, i)
            znum += _num(b, i)
        else
            z += _set(b, i)
        end
    end
    z += znum
    z = _cut(z, _set(b, k), b.v, b.ϵ_sg, _sparsity(g, k), b.cut, false)
    _store_set!(b, z, k)
end

function fprop_n!(t::Relax, v::Val{MIN}, g::DAT, b::RelaxCache{V,N,T}, k::Int) where {V,N,T<:RelaxTag}
    z = Inf*one(MC{N,T})
    znum = Inf
    for i in _children(g, k)
        if _is_num(b, i)
            znum = min(znum, _num(b, i))
        else
            z = min(z, _set(b, i))
        end
    end
    z = min(z, znum)
    z = _cut(z, _set(b, k), b.v, b.ϵ_sg, _sparsity(g, k), b.cut, false)
    _store_set!(b, z, k)
end

function fprop_n!(t::Relax, v::Val{MAX}, g::DAT, b::RelaxCache{V,N,T}, k::Int) where {V,N,T<:RelaxTag}
    z = -Inf*one(MC{N,T})
    znum = -Inf
    for i in _children(g, k)
        if _is_num(b, i)
            znum = max(znum, _num(b, i))
        else
            z = max(z, _set(b, i))
        end
    end
    z = max(z, znum)
    z = _cut(z, _set(b, k), b.v, b.ϵ_sg, _sparsity(g, k), b.cut, false)
    _store_set!(b, z, k)
end

function fprop_n!(t::Relax, ::Val{MULT}, g::DAT, b::RelaxCache{V,N,T}, k::Int) where {V,N,T<:RelaxTag}
    #println("MULT N")
    dp = b.dp
    dP = b.dP
    s = _sparsity(g, 1)
    #l = b.#Float64[_lbd(b.v, i) for i in s]
    #u = #Float64[_ubd(b.v, i) for i in s]
    #P = Interval.(l, u)
    #p0 = b.v.x0[s]
    #p = b.v.x[s]
    z = one(MC{N,T})
    zi = one(MC{N,T}) 
    znum = one(Float64)
    first_set = true
    for i in _children(g, k)
        if _is_num(b, i)
            znum = znum*_num(b, i)
        else
            x = _set(b, i)
            xi = _info(b, i)
            if b.use_apriori_mul
                if !first_set
                    xr = xi
                    zr = zi

                    xr_cv = cv(xr)
                    zr_cv = cv(zr)
                    xr_cv_grad = cv_grad(xr)
                    zr_cv_grad = cv_grad(zr)
                    
                    xrn_cc = -cc(xr)
                    zrn_cc = -cc(zr)
                    xrn_cc_grad = -cc_grad(xr)
                    zrn_cc_grad = -cc_grad(zr)

                    t3 = affine_expand_del(dP, xr_cv, xr_cv_grad, s)
                    t4 = affine_expand_del(dP, zr_cv, zr_cv_grad, s)
                    s3 = affine_expand_del(dP, xrn_cc, xrn_cc_grad, s)
                    s4 = affine_expand_del(dP, zrn_cc, zrn_cc_grad, s)

                    z = x*z
                    wIntv = z.Intv
                    if (t3 < x.Intv.hi) || (t4 < z.Intv.hi)
                        t1 = affine_expand_del(dp, xr_cv, xr_cv_grad, s)
                        t2 = affine_expand_del(dp, zr_cv, zr_cv_grad, s)
                        za_l = McCormick.mult_apriori_kernel(x, z, wIntv, t1, t2, t3, t4, xr_cv_grad, zr_cv_grad)
                        z = z ∩ za_l
                    end
                    if (s3 < -x.Intv.lo) || (s4 < -z.Intv.lo)
                        s1 = affine_expand_del(dp, xrn_cc, xrn_cc_grad, s)
                        s2 = affine_expand_del(dp, zrn_cc, zrn_cc_grad, s)
                        za_u = McCormick.mult_apriori_kernel(-x, -z, wIntv, s1, s2, s3, s4, xrn_cc_grad, zrn_cc_grad)
                        z = z ∩ za_u
                    end
                    zi = xi*zi
                else
                    first_set = false
                    zi = xi
                    z = x
                end
            else
                z = z*x
            end
        end
    end
    z = z*znum
    z = _cut(z, _set(b, k), b.v, b.ϵ_sg, _sparsity(g, k), b.cut, false)
    _store_set!(b, z, k)
end

for F in (PLUS, MULT, MIN, MAX)
    eval(quote
        function fprop!(t::Relax, v::Val{$F}, g::DAT, b::RelaxCache{V,N,T}, k::Int) where {V,N,T<:RelaxTag}
            if _arity(g, k) == 2
                return fprop_2!(Relax(), Val($F), g, b, k)
            end
            fprop_n!(Relax(), Val($F), g, b, k)
        end
    end)
end
function fprop!(t::Relax, v::Val{POW}, g::DAT, b::RelaxCache{V,N,T}, k::Int) where {V,N,T<:RelaxTag}
    x = _child(g, 1, k)
    y = _child(g, 2, k)
    x_is_num = _is_num(b, x)
    y_is_num = _is_num(b, y)
    if y_is_num && isone(_num(b, y))
        z = _set(b, x)
        _store_set!(b, z, k)
    elseif y_is_num && iszero(_num(b, y))
        _store_set!(b, one(_set(b, x)), k)
    elseif x_is_num && y_is_num
        _store_num!(b, _num(b, x)^_num(b, y), k)
    else
        if !x_is_num && y_is_num
            z = _set(b, x)^_num(b, y)
        elseif x_is_num && !y_is_num
            z = _num(b, x)^_set(b, y)
        elseif !x_is_num && !y_is_num
            z = _set(b, x)^_set(b, y)
        end
        z = _cut(z, _set(b, k), b.v, zero(Float64), _sparsity(g,k), b.cut, b.post)
        _store_set!(b, z, k)
    end
end
function fprop!(t::Relax, v::Val{USER}, g::DAT, b::RelaxCache{V,N,T}, k::Int) where {V,N,T<:RelaxTag}
    f = _user_univariate_operator(g, _first_index(g, k))
    x = _set(b, _child(g, 1, k))
    z = f(x)
    z = _cut(z, _set(b, k), b.v, zero(Float64), _sparsity(g, k), b.cut, b.post)
    _store_set!(b, z, k)
end
function fprop!(t::Relax, v::Val{USERN}, g::DAT, b::RelaxCache{V,N,T}, k::Int) where {V,N,T<:RelaxTag}
    mv = _user_multivariate_operator(g, _first_index(g, k))
    n = _arity(g, k)
    set_input = _set_input(b, n)
    i = 1
    for c in _children(g, k)
        if _is_num(b, c)
            x = _num(b, c)
            if !isinf(x)
                set_input[i] = MC{N,T}(x)
            end
        else
            set_input[i] = _set(b, c)
        end
        i += 1
    end
    z = MOI.eval_objective(mv, set_input)::MC{N,T}
    z = _cut(z, _set(b, k), b.v, zero(Float64), _sparsity(g,k), b.cut, b.post)
    _store_set!(b, z, k)
end

for ft in UNIVARIATE_ATOM_TYPES
    f = UNIVARIATE_ATOM_DICT[ft]
    if f == :user || f == :+ || f == :-
        continue
    end
    eval(quote
        function fprop!(t::Relax, v::Val{$ft}, g::DAT, b::RelaxCache{V,N,T}, k::Int) where {V,N,T<:RelaxTag}
            x = _set(b, _child(g, 1, k))
            z = ($f)(x)
            z = _cut(z, _set(b, k), b.v, zero(Float64), _sparsity(g,k), b.cut, b.post)
            _store_set!(b, z, k)
        end
    end)
end

function fprop!(t::Relax, v::Val{ARH}, g::DAT, b::RelaxCache{V,N,T}, k::Int) where {V,N,T<:RelaxTag}
    x = _child(g, 1, k)
    y = _child(g, 2, k)
    x_is_num = _is_num(b, x)
    y_is_num = _is_num(b, y)
    if y_is_num && isone(_num(b, y))
        z = _set(b, x)
        _store_set!(b, z, k)
    elseif y_is_num && iszero(_num(b, y))
        _store_set!(b, zero(V), k)
    else
        if !x_is_num && y_is_num
            z = arh(_set(b, x), _num(b, y))
        elseif x_is_num && !y_is_num
            z = arh(_num(b, x), _set(b, y))
        elseif !x_is_num && !y_is_num
            z = arh(_set(b, x), _set(b, y))
        end
        z = _cut(z, _set(b, k), b.v, zero(Float64), _sparsity(g,k), b.cut, b.post)
        _store_set!(b, z, k)
    end
end

function fprop!(t::Relax, v::Val{LOWER_BND}, g::DAT, b::RelaxCache{V,N,T}, k::Int) where {V,N,T<:RelaxTag}
    z = _set(b, _child(g, 1, k))
    y = _child(g, 2, k)
    if _is_num(b, y)
        z = lower_bnd(z, _num(b, y))
    end
    z = _cut(z, _set(b, k), b.v, zero(Float64), _sparsity(g,k), b.cut, b.post)
    _store_set!(b, z, k)
end

function fprop!(t::Relax, v::Val{UPPER_BND}, g::DAT, b::RelaxCache{V,N,T}, k::Int) where {V,N,T<:RelaxTag}
    z = _set(b, _child(g, 1, k))
    y = _child(g, 2, k)
    if _is_num(b, y)
        z = upper_bnd(z, _num(b, y))
    end
    z = _cut(z, _set(b, k), b.v, zero(Float64), _sparsity(g,k), b.cut, b.post)
    _store_set!(b, z, k)
end

function fprop!(t::Relax, v::Val{BND}, g::DAT, b::RelaxCache{V,N,T}, k::Int) where {V,N,T<:RelaxTag}
    z = _set(b, _child(g, 1, k))
    y = _child(g, 2, k)
    r = _child(g, 3, k)
    if _is_num(b, y) && _is_num(b, r)
        z = bnd(z, _num(b, y), _num(b,r))
    end
    z = _cut(z, _set(b, k), b.v, zero(Float64), _sparsity(g,k), b.cut, b.post)
    _store_set!(b, z, k)
end

function f_init!(t::Relax, g::DAT, b::RelaxCache)
    if b.use_apriori_mul
        s = _sparsity(g, 1)
        v = b.v
        x = v.x
        x0 = v.x0
        for j in s
            b.dp[j] = x[j] - x0[j]
            b.dP[j] = Interval(_lbd(b, j), _ubd(b, j)) - x0[j]
        end
    end
    for k = _node_count(g):-1:1
        if _is_unlocked(b, k)
            c = _node_class(g, k)
            if c == EXPRESSION
                fprop!(t, Expression(), g, b, k)
            elseif c == VARIABLE
                fprop!(t, Variable(), g, b, k)
            elseif c == SUBEXPRESSION
                fprop!(t, Subexpression(), g, b, k)
            end
        end
        b._info[k] = _set(b, k)
    end
    return
end