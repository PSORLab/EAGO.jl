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

function fprop!(t::RelaxCacheAttribute, vt::Variable, g::DAT, b::RelaxCache{V,N,T}, k) where {V,N,T<:RelaxTag}
    i = first_index(g, k)
    x = val(b, i)
    l = lbd(b, i)
    u = ubd(b, i)
    if l == u
        b[k] = x
    else
        #=
        if !isfinite(x)
            if !isfinite(l) && !isfinite(u)
                x = 0.0
            elseif isfinite(u)
                x = u
            elseif isfinite(l)
                x = l
            end
        end
        =#
        z = varset(MC{N,T}, rev_sparsity(g, i, k), x, x, l, u)
        if !first_eval(t, b)
            z = z ∩ interval(b, k)
        end
        #println("[$k] VARIABLE[$i]... = $z \n")
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
            b[k] = cut(z, set(b,k), b.ic.v, b.ϵ_sg, sparsity(g, k), b.cut, false)
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
            #println("[$k] Minus[$x, $y]... = $z")
            b[k] = cut(z, set(b,k), b.ic.v, b.ϵ_sg, sparsity(g, k), b.cut, false)
            #println("[$k] Minus[$x, $y]... = $ztemp after cut \n")
        else
            b[k] = num(b, x) - num(b, y)
        end
    else
        if is_num(b, x)
            b[k] = -num(b, x)
        else
            z = -set(b, x)
            b[k] = cut(-set(b, x), set(b, k), b.ic.v, b.ϵ_sg, sparsity(g, k), b.cut, false)
        end
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
        b[k] = cut(z, set(b, k), b.ic.v, b.ϵ_sg, sparsity(g, k), b.cut, false)
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
        b[k] = cut(z, set(b, k), b.ic.v, b.ϵ_sg, sparsity(g, k), b.cut, false)
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
        b[k] = cut(z, set(b, k), b.ic.v, b.ϵ_sg, sparsity(g, k), b.cut, false)
    end
end

function fprop_2!(t::Relax, v::Val{MULT}, g::DAT, b::RelaxCache{V,N,T}, k::Int) where {V,N,T<:RelaxTag}

    x = child(g, 1, k)
    y = child(g, 2, k)

    if !xy_num(b, x, y)
        if xyset(b, x, y)
            xv = set(b, x)
            yv = set(b, y)
            if b.use_apriori_mul
                dp = b.dp
                dP = b.dP
                p_rel = b.p_rel
                p_diam = b.p_diam
                s = sparsity(g, 1)
                xr = info(b, x)
                yr = info(b, y)
               # @show xr
               # @show yr
                @show xv
                @show yv
                u1max, u2max, v1nmax, v2nmax = estimator_extrema(xr, yr, s, dP)
                @show u1max, u2max, v1nmax, v2nmax
                z = xv*yv
                @show z
                wIntv = z.Intv
                if (u1max < xv.Intv.hi) || (u2max < yv.Intv.hi)
                    u1cv, u2cv, u1cvg, u2cvg = estimator_under(xv, yv, xr, yr, s, dp, dP, p_rel, p_diam)
                    @show u1cv, u2cv, u1cvg, u2cvg
                    za_l = McCormick.mult_apriori_kernel(xv, yv, wIntv, u1cv, u2cv, u1max, u2max, u1cvg, u2cvg)
                    @show za_l
                    z = z ∩ za_l
                end
                if (v1nmax > -xv.Intv.lo) || (v2nmax > -yv.Intv.lo)
                    v1ccn, v2ccn, v1ccgn, v2ccgn = estimator_over(xv, yv, xr, yr, s, dp, dP, p_rel, p_diam)
                    @show v1ccn, v2ccn, v1ccgn, v2ccgn
                    za_u = McCormick.mult_apriori_kernel(-xv, -yv, wIntv, v1ccn, v2ccn, v1nmax, v2nmax, v1ccgn, v2ccgn)
                    @show za_u
                    z = z ∩ za_u
                end
            else
                z = xv*yv
            end
        elseif xset_ynum(b, x, y)
            z = set(b, x)*num(b, y)
        else
            z = num(b, x)*set(b, y)
        end
        #println("[$k] MULT[$x, $y]... = $z (vs. $(b._info[k])) \n")
        b[k] = cut(z, set(b,k), b.ic.v, b.ϵ_sg, sparsity(g, k), b.cut, false)
    else
        b[k] = num(b, x)*num(b, y)
    end
end

function fprop_n!(t::Relax, ::Val{MULT}, g::DAT, b::RelaxCache{V,N,T}, k::Int) where {V,N,T<:RelaxTag}

    z = one(MC{N,T})
    znum = one(Float64)
    numval = true
    if b.use_apriori_mul
        zr = one(V)
        dp = b.dp
        dP = b.dP
        s = sparsity(g, 1)
        for (q,i) in enumerate(children(g, k))
            if is_num(b, i)
                znum = znum*num(b, i)
            else
                numval = false
                x = set(b, i)
                xr = info(b, i)
                u1max, u2max, v1nmax, v2nmax = estimator_extrema(zr, xr, s, dP)
                zv = z*x
                wIntv = zv.Intv
                if (u1max < z.Intv.hi) || (u2max < x.Intv.hi)
                    u1cv, u2cv, u1cvg, u2cvg = estimator_under(zr, xr, s, dp, dP)
                    za_l = McCormick.mult_apriori_kernel(z, x, wIntv, u1cv, u2cv, u1max, u2max, u1cvg, u2cvg)
                    zv = zv ∩ za_l
                end
                if (v1nmax > -z.Intv.lo) || (v2nmax > -x.Intv.lo)
                    v1ccn, v2ccn, v1ccgn, v2ccgn = estimator_under(zr, xr, s, dp, dP)
                    za_u = McCormick.mult_apriori_kernel(-z, -x, wIntv, v1ccn, v2ccn, v1nmax, v2nmax, v1ccgn, v2ccgn)
                    zv = zv ∩ za_u
                end
                zr = zr*xr
                zv = cut(zv, zv, b.ic.v, b.ϵ_sg, sparsity(g, i), b.cut, false)
                z = zv
            end
        end
    else
        for (q,i) in enumerate(children(g, k))
            if is_num(b, i)
                znum = znum*num(b, i)
            else
                numval = false
                x = set(b, i)
                z = z*x
                z = cut(z, z, b.ic.v, b.ϵ_sg, sparsity(g, i), b.cut, false)
            end
        end
    end
    if numval
        b[k] = znum
    else
        z = z*znum
        b[k] = z
    end
end

for F in (PLUS, MULT, MIN, MAX)
    @eval function fprop!(t::Relax, v::Val{$F}, g::DAT, b::RelaxCache{V,N,T}, k) where {V,N,T<:RelaxTag}
        is_binary(g, k) ? fprop_2!(Relax(), Val($F), g, b, k) : fprop_n!(Relax(), Val($F), g, b, k)
    end
end
function fprop!(t::Relax, v::Val{DIV}, g::DAT, b::RelaxCache{V,N,T}, k) where {V,N,T<:RelaxTag}
    fprop_2!(Relax(), v, g, b, k)
end

function fprop!(t::Relax, v::Val{POW}, g::DAT, b::RelaxCache{V,N,T}, k::Int) where {V,N,T<:RelaxTag}
    x = child(g, 1, k)
    y = child(g, 2, k)
    if is_num(b, y) && isone(num(b, y))
        b[k] = set(b, x)
    elseif is_num(b,y) && iszero(num(b, y))
        b[k] = one(MC{N,T})
    elseif !xy_num(b, x, y)
        if xyset(b, x, y)
            z = set(b, x)^set(b, y)
        elseif xset_ynum(b, x, y)
            z = set(b, x)^num(b, y)
        else
            z = num(b, x)^set(b, y)
        end
        b[k] = cut(z, set(b,k), b.ic.v, b.ϵ_sg, sparsity(g, k), b.cut, false)
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
        b[k] = cut(z, set(b, k), b.ic.v, zero(Float64), sparsity(g, k), b.cut, b.post)
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
        b[k] = cut(z, set(b, k), b.ic.v, zero(Float64), sparsity(g,k), b.cut, b.post)
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
            b[k] = cut(z, set(b, k), b.ic.v, zero(Float64), sparsity(g,k), b.cut, b.post)
        end
    end
end

for (F, f) in ((LOWER_BND, :lower_bnd), (UPPER_BND, :upper_bnd))
    @eval function fprop!(t::Relax, v::Val{$F}, g::DAT, b::RelaxCache{V,N,T}, k) where {V,N,T<:RelaxTag}
        y = child(g, 2, k)
        if is_num(b, y)
            z = set(b, child(g, 1, k))
            z = ($f)(z, num(b, y))
            b[k] = cut(z, set(b, k), b.ic.v, zero(Float64), sparsity(g, k), b.cut, b.post)
        end
    end
end

function fprop!(t::Relax, v::Val{BND}, g::DAT, b::RelaxCache{V,N,T}, k) where {V,N,T<:RelaxTag}
    z = set(b, child(g, 1, k))
    y = child(g, 2, k)
    r = child(g, 3, k)
    if is_num(b, y) && is_num(b, r)
        z = bnd(z, num(b, y),num(b, r))
    end
    b[k] = cut(z, set(b, k), b.ic.v, zero(Float64), sparsity(g,k), b.cut, b.post)
end

function f_init!(t::Relax, g::DAT, b::RelaxCache)
    for k = node_count(g):-1:1
        c = node_class(g, k)
        (c == EXPRESSION)    && fprop!(t, Expression(), g, b, k)
        (c == VARIABLE)      && fprop!(t, Variable(), g, b, k)
        (c == SUBEXPRESSION) && fprop!(t, Subexpression(), g, b, k)
        b._info[k] = set(b, k)
    end
    nothing
end