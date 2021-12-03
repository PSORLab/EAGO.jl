
function fprop!(::RelaxInterval, vt::Variable, g::DAT, b::IntervalCache{T}, k) where T<:Real
    i = first_index(g, k)
    l = lbd(b, i)
    u = ubd(b, i)
    b[k] = (l == u) ? Interval(l) : Interval(l, u)
    nothing
end

function fprop!(t::RelaxInterval, ex::Subexpression, g::DAT, b::IntervalCache{T}, k) where T<:Real 
    b[k] = subexpression_set(t, b, first_index(g, k))
end

for (F, f) in ((DIV, :/), (ARH, :arh), (POW, :^))
    @eval fprop!(t::RelaxInterval, v::Val{$F}, g::DAT, b::IntervalCache{T}, k) where T<:Real = (b[k] = ($f)(set(t, b, child(g, 1, k)), set(b, child(g, 2, k))); nothing)
end

function fprop!(t::RelaxInterval, v::Val{MINUS}, g::DAT, b::IntervalCache{T}, k) where T<:Real 
    x = child(g, 1, k)
    b[k] = is_binary(g, k) ? (set(t, b, x) - set(t, b, child(g, 2, k))) : - set(t, b, x)
end

for (F, f) in ((PLUS, :sum), (MIN, :minimum), (MAX, :maximum), (MULT, :prod))
    @eval fprop!(t::Interval, v::Val{$F}, g::DAT, b::IntervalCache{T}, k) where T<:Real = (b[k] = ($f)(i -> set(t, b, i), children(g, k)); nothing)
end

function fprop!(t::RelaxInterval, v::Val{USER}, g::DAT, b::IntervalCache{T}, k) where T<:Real 
    f = user_univariate_operator(g, first_index(g, k))
    b[k] = f(set(t, b, child(g, 1, k)))
end

function fprop!(t::RelaxInterval, v::Val{USERN}, g::DAT, b::IntervalCache{T}, k) where T<:Real
    mv = user_multivariate_operator(g, first_index(g, k))
    set_input = set_input(t, b, arity(g, k))
    for c in children(g, k)
        set_input[i] = set(t, b, c)
    end
    b[k] = MOI.eval_objective(mv, set_input)::Interval{Float64}
end

for ft in UNIVARIATE_ATOM_TYPES
    f = UNIVARIATE_ATOM_DICT[ft]
    (f == :user || f == :+ || f == :-) && continue
    @eval function fprop!(t::RelaxInterval, v::Val{$ft}, g::DAT, b::IntervalCache{T}, k) where T<:Real
        x = child(g, 1, k)
        b[k] = ($f)(set(t, b, x))
    end
end

for (F, f) in ((LOWER_BND, :lower_bnd), (UPPER_BND, :upper_bnd))
    @eval function fprop!(t::RelaxInterval, v::Val{$F}, g::DAT, b::IntervalCache{T}, k) where T<:Real 
        y = child(g, 2, k)
        if is_num(t, b, y)
            z = set(t, b, child(g, 1, k))
            b[k] = ($f)(z, num(t, b, y))
        end
        nothing
    end
end

function fprop!(t::RelaxInterval, v::Val{BND}, g::DAT, b::IntervalCache{T}, k) where T<:Real 
    y = child(g, 2, k)
    r = child(g, 3, k)
    if is_num(t, b, y) && is_num(t, b, r)
        z = set(t, b, child(g, 1, k))
        b[k] = bnd(z, num(t, b, y), num(t, b, r))
    end
    nothing
end