#=
using JuMP
using EAGO

m = Model()

register_eago_operators!(m)
@variable(m, x[1:m])
@auxillary_variable(m, v[1:n])
=#

"""
    ImplicitFunction{H}

Specifies the implicit function y(x) defined by the equality constraint
h(v, x) = 0. Restrictions on the domain of v and x are specified via
constraints provided by JuMP.
```
@constraint(m, (v, x) in ImplicitFunction(h))
```
"""
struct ImplicitFunction{H}
    h::H
    n::Int # Dimension of x
    m::Int # Dimension of v
end

"""
    MIMOExpr{H}

Specifies an expression that defines `v` as an inplace function `h!(v, x)` of `x`.
```
@expression(m, MIMOExpr(h!, v, x))
```
"""
struct MIMOExpr{H}
    h!::H
end
