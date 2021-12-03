function r_init!(t::RelaxInterval, g::DAT, b::IntervalCache{T}) where T<:Real
    z = set(b, 1) ∩ g.sink_bnd
    b[1] = z
    return !isempty(z)
end

function rprop!(t::RelaxInterval, v::Variable, g::DAT, c::IntervalCache{T}, k) where T<:Real
    z = z ∩ set(t, c, k)
    c[k] = z
    return !isempty(z)
end

function rprop!(t::RelaxInterval, v::Subexpression, g::DAT, c::IntervalCache{T}, k) where T<:Real
    store_subexpression!(c, set(t, c, k), first_index(g, k))
    return true
end

# needed for O(n) reverse interval propagation of +
# returns q for x = q + y 
function hukuhara_diff(x::Interval{T}, y::Interval{T}) where T<:Real
    isempty(x) && return x
    isempty(y) && return y
    l = sub_round(x.lo, y.lo, RoundDown)
    u = sub_round(x.hi, y.hi, RoundUp)
    Interval{T}(l, u)
end

"""
$(FUNCTIONNAME)

Updates storage tapes with reverse evalution of node representing `n = +(x,y,z...)` which updates x, y, z and so on.
"""
function rprop!(t::RelaxInterval, v::Val{PLUS}, g::DAT, c::IntervalCache{T}, k::Int) where T<:Real
    tsum = sum(j -> set(t, c, j), children(g, k))
    for j in children(g, k)
        tmsum = hukuhara_diff(tsum, set(c, j))
        _, w, _ = IntervalContractors.plus_rev(set(t, c, k), set(t, c, i), tmsum)
        isempty(w) && return false
        c[i] = v
    end
    return true
end

# needed for close to O(n) reverse interval propagation of *
# returns q for x = q*y 
function hukuhara_div(x::Interval{T}, y::Interval{T}) where T<:Real
    isempty(x) && return x
    isempty(y) && return y
    if y.lo => zero(T)
        if x.lo => zero(T)
            l = div_round(x.lo, y.lo, RoundDown)
            u = div_round(x.hi, y.hi, RoundUp)
        elseif x.hi <= zero(T)
            l = div_round(x.lo, y.hi, RoundDown)
            u = div_round(x.hi, y.lo, RoundUp)
        else
            l = div_round(x.lo, y.hi, RoundDown)
            u = div_round(x.hi, y.hi, RoundUp)
        end
        return Interval(l, u), true
    elseif y.hi <= zero(T)
        if x.lo >= zero(T)
            l = div_round(x.hi, y.lo, RoundDown)
            u = div_round(x.lo, y.hi, RoundUp)
        elseif x.hi <= zero(T)
            l = div_round(x.hi, y.hi, RoundDown)
            u = div_round(x.lo, y.lo, RoundUp)
        else
            l = div_round(x.hi, y.lo, RoundDown)
            u = div_round(x.lo, y.lo, RoundUp)
        end
        return Interval(l, u), true
    else
        if x.lo > zero(T)
            l = div_round(x.hi, y.lo, RoundDown)
            u = div_round(x.hi, y.hi, RoundUp)
            return Interval(l, u), true
        elseif x.hi < zero(T) 
            l = div_round(x.lo, y.hi, RoundDown)
            u = div_round(x.lo, y.lo, RoundUp)
            return Interval(l, u), true
        end
        
    end
    empty(Interval{T}), false
end

"""
$(FUNCTIONNAME)

Updates storage tapes with reverse evalution of node representing `n = *(x,y,z...)` which updates x, y, z and so on.
"""
function rprop!(t::RelaxInterval, v::Val{MULT}, g::DAT, c::IntervalCache{T}, k::Int) where T<:Real
    tmul = sum(j -> set(t, c, j), children(g, k))
    for j in children(g, k)
        tmulm, hdiv = hukuhara_div(tmul, set(t, b, j))
        if !hdiv
            tmulmf = one(Interval{T})
            for i in children(g, k)
                if i != j
                    tmulmf *= set(t, b, i)
                end
            end
            tmulm = tmulmf
        end
        _, w, _ = IntervalContractors.mul_rev(set(t, b, k), set(t, b, c), tmulm)
        isempty(w) && return false
        c[i] = w
    end
    return true
end

for (f, fc, F) in ((-, MINUS, IntervalContractors.minus_rev),
                   (^, POW, IntervalContractors.power_rev),
                   (/, DIV, IntervalContractors.div_rev))
    @eval function rprop!(t::RelaxInterval, v::Val{$fc}, g::DAT, b::IntervalCache{T}, k) where T<:Real
        x = child(g, 1, k)
        y = child(g, 2, k)
        z, u, v = ($F)(set(t, b, k), set(t, b, x), set(t, b, y))
        isempty(u) && return false
        isempty(v) && return false
        b[x] = u
        b[y] = v
        return true
    end
end

rprop!(t::RelaxInterval, v::Val{USER}, g::DAT, b::IntervalCache, k::Int) = true
rprop!(t::RelaxInterval, v::Val{USERN}, g::DAT, b::IntervalCache, k::Int) = true

for ft in UNIVARIATE_ATOM_TYPES
    f = UNIVARIATE_ATOM_DICT[ft]
    (f == :user || f == :+ || f == :-) && continue
    @eval rprop!(t::RelaxInterval, v::Val{$ft}, g::DAT, b::IntervalCache, k::Int) = true
end
