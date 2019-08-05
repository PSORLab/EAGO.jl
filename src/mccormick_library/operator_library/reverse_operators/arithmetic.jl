"""
    plus_rev

Creates reverse McCormick contractor for `a` = `b` +`c`
"""
function plus_rev(a::MC, b::MC, c::MC)  # a = b + c
    bold = b
    b = b ∩ (a - c)
    c = c ∩ (a - bold)
    a, b, c
end
#plus_rev(a,b,c) = plus_rev(promote(a,b,c)...)

"""
    minus_rev

Creates reverse McCormick contractor for `a` = `b`- `c`
"""
function minus_rev(a::MC, b::MC, c::MC)  # a = b - c
    bold = b
    b = b ∩ (a + c)
    c = c ∩ (bold - a)
    a, b, c
end
function minus_rev(a::MC, b::MC)  # a = -b
    b = b ∩ -a
    return a, b
end
#minus_rev(a,b,c) = minus_rev(promote(a,b,c)...)
#minus_rev(a,b) = minus_rev(promote(a,b)...)

"""
    mul_rev

Creates reverse McCormick contractor for `a` = `b`*`c`
"""
function mul_rev(a::MC, b::MC, c::MC)  # a = b * c
    #=
    bflag = (0.0 ∉ b.Intv)
    if bflag
        temp1 = a / b
        ((0.0 ∉ a.Intv) || bflag) && (c = c ∩ temp1)
    end
    cflag = (0.0 ∉ c.Intv)
    if cflag
        temp2 = a / c
        ((0.0 ∉ a.Intv) || cflag) && (b = b ∩ temp2)
    end
    =#
    a,b,c
end
mul_rev(a::MC{N},b::MC{N},c::Float64) where N = mul_rev(a,b,MC{N}(c))
mul_rev(a::MC{N},b::Float64,c::MC{N}) where N = mul_rev(a,MC{N}(b),c)

"""
    div_rev

Creates reverse McCormick contractor for `a` = `b`/`c`
"""
function div_rev(a::MC, b::MC, c::MC)  # a = b / c
    #=
    b = b ∩ (a * c)
    if (~isempty(b) && ~in(0.0, a.Intv))
        c = c ∩ (b / a)
    end
    =#
    a,b,c
end
#div_rev(a,b,c) = div_rev(promote(a,b,c)...)

"""
    inv_rev

Creates reverse McCormick contractor for `a` = `inv(b)`
"""
function inv_rev(a::MC, b::MC)  # a = inv(b)
    ~in(0.0, a.Intv) && (b = b ∩ inv(a))
    a,b
end
inv_rev(a,b) = inv_rev(promote(a,b)...)

"""
    power_rev

Creates reverse McCormick contractor for `a` = `b`^`c`
"""
function power_rev(a::MC, b::MC, c::MC)  # a = b^c
    #=
    if (~isempty(b) && ~isempty(c))
        ~in(0.0, c.Intv) && (b = b ∩ (a^(inv(c))))
        if ~in(0.0, a.Intv)
            if (0.0 < b.Intv.lo < Inf)
                blog = log(b)
                if ~in(0.0, blog.Intv)
                    c = c ∩ (log(a) / blog)
                end
            end
        end
    else
        a = empty(a)
    end
    =#
    a,b,c
end


"""
    sqrt_rev

Creates reverse McCormick contractor for `a` = `sqrt(b)`
"""
function sqrt_rev(a::MC, b::MC)  # a = sqrt(b)
    #b = b ∩ (a^2)
    a,b
end
#sqr_rev(f, x)  = power_rev(f,x,2)


"""
    abs_rev!

Creates reverse McCormick contractor for `a` = `abs(b)`
"""
abs_rev(a::MC, b::MC) = (a,b)
#=
function abs_rev!(y::MC{N}, x::MC{N}) where N   # y = abs(x); refine x

    y_new = y ∩ (0..∞)

    x1 = y_new ∩ x
    x2 = -(y_new ∩ (-x))

    cv =
    cc =
    cv_grad =
    cc_grad =
    Intv = hull(Intv(x1),Intv(x2))


    y = MC{N}(cv, cc, Intv, cv_grad, cc_grad, y.cnst)

    return
end
=#
