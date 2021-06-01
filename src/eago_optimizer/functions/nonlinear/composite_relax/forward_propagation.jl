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

@inline affine_expand(x, x0, fx0, ∇fx0) = fx0 + dot(∇fx0, x - x0)

f_init!(::Relax, g::ALLGRAPHS, b::RelaxCache) = nothing

function _var_set(::Type{MC{N,T}}, i::Int, x_cv::Float64, x_cc::Float64, l::Float64, u::Float64) where {N,T<:RelaxTag}
    v = seed_gradient(i, Val(N))
    return MC{N,T}(x_cv, x_cc, Interval{Float64}(l, u), v, v, false)
end

function fprop!(t::Relax, vt::Variable, g::ALLGRAPHS, b::RelaxCache{N,T}, k::Int) where {N,T<:RelaxTag}
    i = _first_index(g, k)
    x = _val(b, i)
    z = _var_set(MC{N,T}, _rev_sparsity(g, i, k), x, x, _lbd(b, i), _ubd(b, i))
    if !_first_eval(b)
        z = z ∩ _interval(b, k)
    end
    _store_set!(b, z, k)
    (b.first_eval && b.use_apriori_mul) && _store_info!(b, z, k)
    return
end

function expand_set(::Type{MC{N2,T}}, x::MC{N1,T}, fsparse::Vector{Int},
                    subsparse::Vector{Int}, cv_buffer::Vector{Float64},
                    cc_buffer::Vector{Float64}) where {N1, N2, T<:RelaxTag}

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

function fprop!(t::Relax, ex::Subexpression, g::ALLGRAPHS, c::RelaxCache{N,T}, k::Int) where {N,T<:RelaxTag}
    d = _subexpression_value(c, _first_index(g, k))
    z = expand_set(MC{N,T}, d.set[1], _sparsity(g, k), _sparsity(sub), c.cv_buffer, c.cc_buffer)
    _store_set!(c, z, k)
    (b.first_eval && b.use_apriori_mul) && _store_info!(b, z, k)
    return
end

function set_value_post(z::MC{N,T}, v::VariableValues{Float64}, s::Vector{Int}, ϵ::Float64) where {N,T<:RelaxTag}
    lower = z.cv
    upper = z.cc
    lower_refinement = true
    upper_refinement = true
    @inbounds for i = 1:N
        cv_val = z.cv_grad[i]
        cc_val = z.cc_grad[i]
        i_sol = s[i]
        x_z = v.x[i_sol]
        lower_bound = v.lbd[i_sol]
        upper_bound = v.ubd[i_sol]
        if lower_refinement
            if cv_val > zero(Float64)
                if isinf(lower_bound)
                    !upper_refinement && break
                    lower_refinement = false
                else
                    lower += cv_val*(lower_bound - x_z)
                end
            else
                if isinf(upper_bound)
                    !upper_refinement && break
                    lower_refinement = false
                else
                    lower += cv_val*(upper_bound - x_z)
                end
            end
        end
        if upper_refinement
            if cc_val > zero(Float64)
                if isinf(lower_bound)
                    !lower_refinement && break
                    upper_refinement = false
                else
                    upper += cc_val*(upper_bound - x_z)
                end
            else
                if isinf(upper_bound)
                    !lower_refinement && break
                    upper_refinement = false
                else
                    upper += cc_val*(lower_bound - x_z)
                end
            end
        end
    end

    if lower_refinement && (z.Intv.lo + ϵ > lower)
        lower = z.Intv.lo
    elseif !lower_refinement
        lower = z.Intv.lo
    else
        lower -= ϵ
    end

    if upper_refinement && (z.Intv.hi - ϵ < upper)
        upper = z.Intv.hi
    elseif !upper_refinement
        upper = z.Intv.hi
    else
        upper += ϵ
    end

    return MC{N,T}(z.cv, z.cc, Interval{Float64}(lower, upper), z.cv_grad, z.cc_grad, z.cnst)
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
function _cut(x::MC{N,T}, lastx::MC{N,T}, v::VariableValues, ϵ::Float64, s::Vector{Int},
              post::Bool, cut::Bool, cut_interval::Bool) where {N,T<:RelaxTag}

    if post && cut && cut_interval
        return set_value_post(x ∩ lastx.Intv, v, s, ϵ)
    elseif post && !cut
        return set_value_post(x, v, s, ϵ)
    elseif !post && cut && cut_interval
        return x ∩ lastx.Intv
    elseif !post && cut && !cut_interval
        return x ∩ lastx
    end
    return x
end

for (f, F, fc) in ((:fprop_2!, PLUS, :+), (:fprop_2!, MIN, :min), (:fprop_2!, MAX, :max), (:fprop!, DIV, :/))
    eval(quote
        function ($f)(t::Relax, v::Val{$F}, g::ALLGRAPHS, b::RelaxCache{N,T}, k::Int) where {N,T<:RelaxTag}
            x = _child(g, 1, k)
            y = _child(g, 2, k)
            if !_is_num(b, x) && _is_num(b, y)
                z = ($fc)(_set(b, x), _num(b, y))
            elseif _is_num(b, x) && !_is_num(b, y)
                z = ($fc)(_num(b, x), _set(b, y))
            else
                z = ($fc)(_set(b, x), _set(b, y))
            end
            z = _cut(z, _set(b,k), b.v, b.ϵ_sg, _sparsity(g, k), false, b.cut, b.cut_interval)
            _store_set!(b, z, k)
            (b.first_eval && b.use_apriori_mul) && _store_info!(b, z, k)
            return
        end
    end)
end
function fprop!(t::Relax, v::Val{MINUS}, g::ALLGRAPHS, b::RelaxCache{N,T}, k::Int) where {N,T<:RelaxTag}
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
    z = _cut(z, _set(b, k), b.v, b.ϵ_sg, _sparsity(g, k), false, b.cut, b.cut_interval)
    _store_set!(b, z, k)
    (b.first_eval && b.use_apriori_mul) && _store_info!(b, z, k)
    return
end


function fprop_2!(t::Relax, v::Val{MULT}, g::ALLGRAPHS, b::RelaxCache{N,T}, k::Int) where {N,T<:RelaxTag}
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
        if b.use_apriori_mul
            xr = _info(b, x)
            yr = _info(b, y)
            z = mult_apriori_kernel(xv, yv, xv.Intv*yv.Intv,
                                    affine_expand(p, p0, xr.cv, xr.cv_grad),
                                    affine_expand(p, p0, yr.cv, yr.cv_grad),
                                    hi(affine_expand(P, p0, xr.cv, xr.cv_grad)),
                                    hi(affine_expand(P, p0, yr.cv, yr.cv_grad)),
                                    affine_expand(p, p0, xr.cc, xr.cc_grad),
                                    affine_expand(p, p0, yr.cc, yr.cc_grad),
                                    lo(affine_expand(P, p0, xr.cc, xr.cc_grad)),
                                    lo(affine_expand(P, p0, yr.cc, yr.cc_grad)),
                                    xr.cv_grad, yr.cv_grad, xr.cc_grad, yr.cc_grad)
        else
            z = xv*yv
        end
    end
    z = _cut(z, _set(b,k), b.v, b.ϵ_sg, _sparsity(g, k), false, b.cut, b.cut_interval)
    _store_set!(b, z, k)
    (b.first_eval && b.use_apriori_mul) && _store_info!(b, z, k)
    return
end

for (F, FT, SV, NV) in ((:+,   PLUS, :(zero(MC{N,T})), zero(Float64)),
                        (:min, MIN,  :(inf(MC{N,T})),  Inf),
                        (:max, MAX,  :(-inf(MC{N,T})), -Inf))
    eval(quote 
            function fprop_n!(t::Relax, v::Val{$FT}, g::ALLGRAPHS, b::RelaxCache{N,T}, k::Int) where {N,T<:RelaxTag}
                z = $SV
                znum = $NV
                for i in _children(g, k)
                    if _is_num(b, i)
                        znum = ($F)(znum, _num(b, i))
                        continue
                    end
                    z = ($F)(z, _set(b, i))
                end
                z = ($F)(z, znum)
                z = _cut(z, _set(b, k), b.v, b.ϵ_sg, _sparsity(g, k), false, b.cut, b.cut_interval)
                _store_set!(b, z, k)
                (b.first_eval && b.use_apriori_mul) && _store_info!(b, z, k)
                return
            end
        end)
end
function fprop_n!(t::Relax, ::Val{MULT}, g::ALLGRAPHS, b::RelaxCache{N,T}, k::Int) where {N,T<:RelaxTag}
    z = one(MC{N,T})
    znum = one(Float64)
    count = 0
    for i in _children(g, k)
        if _is_num(b, i)
            znum = znum*_num(b, i)
    #        continue
    #    end
    #=
        if b.use_apriori_mul
            if count == 0
                zaff = zaff*_info(b,i)::MC{N,T}
                z = z*_set(b, i)::MC{N,T}
            else
                xi = _set(b, i)::MC{N,T}
                xref = _info(b,i)::MC{N,T}
                z = mult_apriori_kernel(z, xi, z.Intv*xi.Intv,
                                        affine_expand(p, p0, zaff.cv, zaff.cv_grad),
                                        affine_expand(p, p0, xref.cv, xref.cv_grad),
                                        hi(affine_expand(P, p0, zaff.cv, zaff.cv_grad)),
                                        hi(affine_expand(P, p0, xref.cv, xref.cv_grad)),
                                        affine_expand(p, p0, zaff.cc, zaff.cc_grad),
                                        affine_expand(p, p0, xref.cc, xref.cc_grad),
                                        lo(affine_expand(P, p0, zaff.cc, zaff.cc_grad)),
                                        lo(affine_expand(P, p0, xref.cc, xref.cc_grad)),
                                        zaff.cv_grad, xref.cv_grad, zaff.cc_grad, xref.cc_grad)
                zaff = zaff*xref
            end
            =#
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

for F in (PLUS, MULT, MIN, MAX)
    eval(quote
            function fprop!(t::Relax, v::Val{$F}, g::ALLGRAPHS, b::RelaxCache{N,T}, k::Int) where {N,T<:RelaxTag}
            n = _arity(g, k)
            if n == 2
                return fprop_2!(Relax(), Val($F), g, b, k)
            end
            fprop_n!(Relax(), Val($F), g, b, k)
            end
    end)
end
function fprop!(t::Relax, v::Val{POW}, g::ALLGRAPHS, b::RelaxCache{N,T}, k::Int) where {N,T<:RelaxTag}
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
            z = _set(b, x)^_num(b, y)
        elseif x_is_num && !y_is_num
            z = _num(b, x)^_set(b, y)
        elseif !x_is_num && !y_is_num
            z = _set(b, x)^_set(b, y)
        end
        z = _cut(z, _set(b, k), b.v, zero(Float64), _sparsity(g,k), b.post, b.cut, b.cut_interval)
        _store_set!(b, z, k)
    end
    (b.first_eval && b.use_apriori_mul) && _store_info!(b, z, k)
    return
end
function fprop!(t::Relax, v::Val{USER}, g::ALLGRAPHS, b::RelaxCache{N,T}, k::Int) where {N,T<:RelaxTag}
    f = _user_univariate_operator(g, _index(g, k))
    x = _set(b, _child(g, 1, k))
    z = _cut(f(x), _set(b, k), b.v, zero(S), _sparsity(g, k), b.post, b.cut, b.cut_interval)
    _store_set!(b, z, k)
    (b.first_eval && b.use_apriori_mul) && _store_info!(b, z, k)
    return
end
function fprop!(t::Relax, v::Val{USERN}, g::ALLGRAPHS, b::RelaxCache{N,T}, k::Int) where {N,T<:RelaxTag}
    mv = _user_multivariate_operator(g, _index(g, k))
    n = _arity(g, k)
    set_input = _set_input(b, n)
    i = 1
    for c in _children(g, k)
        if _is_num(b, c)
            x = _num(b, c)
            if !isinf(x)
                @inbounds set_input[i] = MC{N,T}(x)
            end
        else
            @inbounds set_input[i] = _set(b, c)
        end
        i += 1
    end
    z = MOI.eval_objective(mv, set_input)
    z = _cut(z, _set(b, k), b.v, zero(Float64), _sparsity(g,k), b.post, b.cut, b.cut_interval)
    _store_set!(b, z, k)
    (b.first_eval && b.use_apriori_mul) && _store_info!(b, z, k)
    return
end

for ft in UNIVARIATE_ATOM_TYPES
    f = UNIVARIATE_ATOM_DICT[ft]
    if f == :user || f == :+ || f == :-
        continue
    end
    eval(quote
        function fprop!(t::Relax, v::Val{$ft}, g::ALLGRAPHS, b::RelaxCache{N,T}, k::Int) where {N,T<:RelaxTag}
            x = _set(b, _child(g, 1, k))
            z = ($f)(x)
            z = _cut(z, _set(b, k), b.v, zero(Float64), _sparsity(g,k), b.post, b.cut, b.cut_interval)
            _store_set!(b, z, k)
            (b.first_eval && b.use_apriori_mul) && _store_info!(b, z, k)
            return
        end
    end)
end

function fprop!(t::Relax, v::Val{ARH}, g::ALLGRAPHS, b::RelaxCache{N,T}, k::Int) where {N,T<:RelaxTag}
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
        z = _cut(z, _set(b, k), b.v, zero(Float64), _sparsity(g,k), b.post, b.cut, b.cut_interval)
        _store_set!(b, z, k)
    end
    (b.first_eval && b.use_apriori_mul) && _store_info!(b, z, k)
    return
end

function fprop!(t::Relax, v::Val{LOWER_BND}, g::ALLGRAPHS, b::RelaxCache{N,T}, k::Int) where {N,T<:RelaxTag}
    z = _set(b, _child(g, 1, k))
    y = _child(g, 2, k)
    if _is_num(b, y)
        z = lower_bnd(z, _num(b, y))
    end
    z = _cut(z, _set(b, k), b.v, zero(Float64), _sparsity(g,k), b.post, b.cut, b.cut_interval)
    _store_set!(b, z, k)
    (b.first_eval && b.use_apriori_mul) && _store_info!(b, z, k)
    return
end

function fprop!(t::Relax, v::Val{UPPER_BND}, g::ALLGRAPHS, b::RelaxCache{N,T}, k::Int) where {N,T<:RelaxTag}
    z = _set(b, _child(g, 1, k))
    y = _child(g, 2, k)
    if _is_num(b, y)
        z = upper_bnd(z, _num(b, y))
    end
    z = _cut(z, _set(b, k), b.v, zero(Float64), _sparsity(g,k), b.post, b.cut, b.cut_interval)
    _store_set!(b, z, k)
    (b.first_eval && b.use_apriori_mul) && _store_info!(b, z, k)
    return
end

function fprop!(t::Relax, v::Val{BND}, g::ALLGRAPHS, b::RelaxCache{N,T}, k::Int) where {N,T<:RelaxTag}
    z = _set(b, _child(g, 1, k))
    y = _child(g, 2, k)
    r = _child(g, 3, k)
    if _is_num(b, y) && _is_num(b, r)
        z = bnd(z, _num(b, y), _num(b,r))
    end
    z = _cut(z, _set(b, k), b.v, zero(Float64), _sparsity(g,k), b.post, b.cut, b.cut_interval)
    _store_set!(b, z, k)
    (b.first_eval && b.use_apriori_mul) && _store_info!(b, z, k)
    return
end
