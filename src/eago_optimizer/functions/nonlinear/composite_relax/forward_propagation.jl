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

xnum_yset(b, x, y) = is_num(b, x) && !is_num(b, y)
xset_ynum(b, x, y) = !is_num(b, x) && is_num(b, y)
xy_num(b, x, y) = is_num(b, x) && is_num(b, y)
xyset(b, x, y) = !(is_num(b, x) || is_num(b, y))

function varset(::Type{MC{N,T}}, i, x_cv, x_cc, l, u) where {V,N,T<:RelaxTag}
    v = seed_gradient(i, Val(N))
    return MC{N,T}(x_cv, x_cc, Interval{Float64}(l, u), v, v, false)
end

function fprop!(t::Relax, vt::Variable, g::DAT, b::RelaxCache{V,N,T}, k) where {V,N,T<:RelaxTag}
    i = first_index(g, k)
    x = val(b, i)
    l = lbd(b, i)
    u = ubd(b, i)
    if l == u
        b[k] = x
    else
        z = varset(MC{N,T}, rev_sparsity(g, i, k), x, x, l, u)
        if !first_eval(t, b)
            z = z ∩ interval(b, k)
        end
        b[k] = z
    end
    nothing
end

function fprop!(t::Relax, ex::Subexpression, g::DAT, b::RelaxCache{V,N,T}, k) where {V,N,T<:RelaxTag}
    x =  first_index(g, k)
    if subexpression_is_num(b, x)
        b[k] = subexpression_num(b, x)
    else
        b[k] = subexpression_set(b, x)
    end
end

for (F, f) in ((PLUS, :+), (MIN, :min), (MAX, :max), (DIV, :/), (ARH, :arh))
    @eval function fprop_2!(t::Relax, v::Val{$F}, g::DAT, b::RelaxCache{V,N,T}, k) where {V,N,T<:RelaxTag}
        x = child(g, 1, k)
        y = child(g, 2, k)
        if !xy_num(b, x, y)
            if xyset(b, x, y)
                z = ($f)(set(b, x), set(b, y))
            elseif xset_ynum(b, x, y)
                z = ($f)(set(b, x), num(b, y))
            else
                z = ($f)(num(b, x), set(b, y))
            end
            b[k] = cut(z, set(b,k), b.v, b.ϵ_sg, sparsity(g, k), b.cut, false)
        else
            b[k] = ($f)(num(b, x), num(b, y))
        end
    end
end

function fprop!(t::Relax, v::Val{MINUS}, g::DAT, b::RelaxCache{V,N,T}, k) where {V,N,T<:RelaxTag}
    x = child(g, 1, k)
    if is_binary(g, k)
        y = child(g, 2, k)
        if !xy_num(b, x, y)
            if xyset(b, x, y)
                z = set(b, x) - set(b, y)
            elseif xset_ynum(b, x, y)
                z = set(b, x) - num(b, y)
            else
                z = num(b, x) - set(b, y)
            end
            b[k] = cut(z, set(b,k), b.v, b.ϵ_sg, sparsity(g, k), b.cut, false)
        else
            b[k] = num(b, x) - num(b, y)
        end
    else
        if is_num(b, x)
            b[k] = -num(b, x)
        else
            z = -set(b, x)
            b[k] = cut(-set(b, x), set(b, k), b.v, b.ϵ_sg, sparsity(g, k), b.cut, false)
        end
    end
end

function fprop_2!(t::Relax, v::Val{MULT}, g::DAT, b::RelaxCache{V,N,T}, k::Int) where {V,N,T<:RelaxTag}

    x = child(g, 1, k)
    y = child(g, 2, k)

    if !xy_num(b, x, y)
        if xyset(b, x, y)
            println("BINARY MULTIPLICATION")
            xv = set(b, x)
            yv = set(b, y)
            if b.use_apriori_mul
                dp = b.dp
                dP = b.dP
                s = sparsity(g, 1)

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
                println("start z = $(z)")
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
                println("finish z = $(z)")
            else
                z = xv*yv
            end
        elseif xset_ynum(b, x, y)
            z = set(b, x)*num(b, y)
        else
            z = num(b, x)*set(b, y)
        end
        b[k] = cut(z, set(b,k), b.v, b.ϵ_sg, sparsity(g, k), b.cut, false)
    else
        b[k] = num(b, x)*num(b, y)
    end
end

function fprop_n!(t::Relax, v::Val{PLUS}, g::DAT, b::RelaxCache{V,N,T}, k::Int) where {V,N,T<:RelaxTag}
    z = zero(MC{N,T})
    znum = 0.0
    numval = true
    for i in children(g, k)
        if is_num(b, i)
            znum += num(b, i)
        else
            numval = false
            z += set(b, i)
        end
    end
    if numval
        b[k] = znum
    else
        z += znum
        b[k] = cut(z, set(b, k), b.v, b.ϵ_sg, sparsity(g, k), b.cut, false)
    end
end

function fprop_n!(t::Relax, v::Val{MIN}, g::DAT, b::RelaxCache{V,N,T}, k::Int) where {V,N,T<:RelaxTag}
    z = Inf*one(MC{N,T})
    znum = Inf
    numval = true
    for i in children(g, k)
        if is_num(b, i)
            znum = min(znum, num(b, i))
        else
            numval = false
            z = min(z, set(b, i))
        end
    end
    if numval
        b[k] = znum
    else
        z = min(z, znum)
        b[k] = cut(z, set(b, k), b.v, b.ϵ_sg, sparsity(g, k), b.cut, false)
    end
end

function fprop_n!(t::Relax, v::Val{MAX}, g::DAT, b::RelaxCache{V,N,T}, k::Int) where {V,N,T<:RelaxTag}
    z = -Inf*one(MC{N,T})
    znum = -Inf
    numval = true
    for i in children(g, k)
        if is_num(b, i)
            znum = max(znum, num(b, i))
        else
            numval = false
            z = max(z, set(b, i))
        end
    end
    if numval
        b[k] = znum
    else
        z = max(z, znum)
        b[k] = cut(z, set(b, k), b.v, b.ϵ_sg, sparsity(g, k), b.cut, false)
    end
end

function fprop_n!(t::Relax, ::Val{MULT}, g::DAT, b::RelaxCache{V,N,T}, k::Int) where {V,N,T<:RelaxTag}
    println(" ")
    println("MUL NARY")
    dp = b.dp
    dP = b.dP
    s = sparsity(g, 1)
    z = one(MC{N,T})
    zr = one(MC{N,T})
    zok = one(MC{N,T})
    znum = one(Float64)
    firstset = true
    numval = true
    for i in children(g, k)
        @show i
        if is_num(b, i)
            znum = znum*num(b, i)
        else
            numval = false
            x = set(b, i)
            xr = info(b, i)
            if b.use_apriori_mul
                if !firstset
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

                    wIntv = x.Intv*z.Intv
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
                    zr = xr*zr
                    zok = x*zok
                else
                    firstset = false
                    zr = xr
                    z = x
                    zok = x
                end
                @show numval
                @show z
                @show zok
                @show znum
                println(" ")
            else
                z = z*x
            end
        end
    end
    if numval
        b[k] = znum
    else
        z = z*znum
        b[k] = cut(z, set(b, k), b.v, b.ϵ_sg, sparsity(g, k), b.cut, false)
    end
end

for F in (PLUS, MULT, MIN, MAX)
    @eval function fprop!(t::Relax, v::Val{$F}, g::DAT, b::RelaxCache{V,N,T}, k) where {V,N,T<:RelaxTag}
        is_binary(g, k) ? fprop_2!(Relax(), Val($F), g, b, k) : fprop_n!(Relax(), Val($F), g, b, k)
    end
end

function fprop!(t::Relax, v::Val{POW}, g::DAT, b::RelaxCache{V,N,T}, k::Int) where {V,N,T<:RelaxTag}
    x = child(g, 1, k)
    y = child(g, 2, k)
    if is_num(b, y) && isone(num(b, y))
        b[k] = set(b, x)
    elseif is_num(b,y) && iszero(_num(b, y))
        b[k] = one(V)
    elseif !xy_num(b, x, y)
        if xyset(b, x, y)
            z = set(b, x)^set(b, y)
        elseif xset_ynum(b, x, y)
            z = set(b, x)^num(b, y)
        else
            z = num(b, x)^set(b, y)
        end
        b[k] = cut(z, set(b,k), b.v, b.ϵ_sg, sparsity(g, k), b.cut, false)
    else
        b[k] = num(b, x)^num(b, y)
    end
end

function fprop!(t::Relax, v::Val{USER}, g::DAT, b::RelaxCache{V,N,T}, k::Int) where {V,N,T<:RelaxTag}
    f = user_univariate_operator(g, first_index(g, k))
    x = child(g, 1, k)
    if is_num(b, x)
        b[k] = f(num(b, x))
    else
        z = f(set(b, x))
        b[k] = cut(z, set(b, k), b.v, zero(Float64), sparsity(g, k), b.cut, b.post)
    end
end

function fprop!(t::Relax, v::Val{USERN}, g::DAT, b::RelaxCache{V,N,T}, k::Int) where {V,N,T<:RelaxTag}
    mv = user_multivariate_operator(g, first_index(g, k))
    n = arity(g, k)
    set_input = _set_input(b, n)
    num_input = _num_input(b, n)
    anysets = false
    i = 1
    for c in children(g, k)
        if is_num(b, c)
            x = num(b, c)
            if !isinf(x)
                set_input[i] = MC{N,T}(x)
                num_input[i] = x
            end
        else
            set_input[i] = set(b, c)
            anysets = true
        end
        i += 1
    end
    if anysets
        z = MOI.eval_objective(mv, set_input)::MC{N,T}
        b[k] = cut(z, set(b, k), b.v, zero(Float64), sparsity(g,k), b.cut, b.post)
    else
        b[k] = MOI.eval_objective(mv, num_input)
    end
end

for ft in UNIVARIATE_ATOM_TYPES
    f = UNIVARIATE_ATOM_DICT[ft]
    (f == :user || f == :+ || f == :-) && continue
    @eval function fprop!(t::Relax, v::Val{$ft}, g::DAT, b::RelaxCache{V,N,T}, k) where {V,N,T<:RelaxTag}
        x = child(g, 1, k)
        if is_num(b, x)
            return b[k] = ($f)(num(b, x))
        else
            z = ($f)(set(b, x))
            b[k] = cut(z, set(b, k), b.v, zero(Float64), sparsity(g,k), b.cut, b.post)
        end
    end
end

for (F, f) in ((LOWER_BND, :lower_bnd), (UPPER_BND, :upper_bnd))
    @eval function fprop!(t::Relax, v::Val{$F}, g::DAT, b::RelaxCache{V,N,T}, k) where {V,N,T<:RelaxTag}
        y = child(g, 2, k)
        if is_num(b, y)
            z = set(b, child(g, 1, k))
            z = ($f)(z, num(b, y))
            b[k] = cut(z, set(b, k), b.v, zero(Float64), sparsity(g, k), b.cut, b.post)
        end
    end
end

function fprop!(t::Relax, v::Val{BND}, g::DAT, b::RelaxCache{V,N,T}, k) where {V,N,T<:RelaxTag}
    z = set(b, child(g, 1, k))
    y = child(g, 2, k)
    r = child(g, 3, k)
    if is_num(b, y) && is_num(b, r)
        z = bnd(z, num(b, y),_num(b, r))
    end
    b[k] = cut(z, set(b, k), b.v, zero(Float64), sparsity(g,k), b.cut, b.post)
end

function f_init!(t::Relax, g::DAT, b::RelaxCache)
    for k = node_count(g):-1:1
        c = node_class(g, k)
        (c == EXPRESSION)    && (fprop!(t, Expression(), g, b, k);    continue)
        (c == VARIABLE)      && (fprop!(t, Variable(), g, b, k);      continue)
        (c == SUBEXPRESSION) && (fprop!(t, Subexpression(), g, b, k); continue)
        b._info[k] = set(b, k)
    end
    nothing
end