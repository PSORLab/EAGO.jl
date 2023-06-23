# Copyright (c) 2018: Matthew Wilhelm, Robert Gottlieb, Dimitri Alston,
# Matthew Stuber, and the University of Connecticut (UConn).
# This code is licensed under the MIT license (see LICENSE.md for full details).
################################################################################
# EAGO
# A development environment for robust and global optimization.
# https://github.com/PSORLab/EAGO.jl
################################################################################
# src/eago_optimizer/unsafe_utilities.jl
# Inbounds and non-allocating versions of simple Julia functions.
################################################################################

"""
$(TYPEDSIGNATURES)

Performs `map!(f, y, x)` in an unsafe manner if y[i] is true, else no-op.
Assumes `n == length(x) == length(y)`. About 2x faster for small arrays (n < 1000).

[Unused]
"""
function unsafe_check_fill!(f, y::Vector{T}, x::T, n::Int) where T
    i = 1
    m = n + 1
    while i < m
        flag = f(@inbounds y[i])
        if flag
            @inbounds y[i] = x
        end
        i += 1
    end
    nothing
end

_rf_findmax((fm, m), (fx, x)) = Base.isless(fm, fx) ? (fx, x) : (fm, m)
map_findmax(f, itr) = mapfoldl(((k, v),) -> (f(v), k), _rf_findmax, pairs(itr))
map_argmax(f, itr) = map_findmax(f, itr)[2]

argmax(f, domain) = mapfoldl(x -> (f(x), x), _rf_findmax, domain)[2]

relative_gap(L::Float64, U::Float64) = ((L > -Inf) && (U < Inf)) ?  abs(U - L)/(max(abs(L), abs(U))) : Inf
relative_tolerance(L::Float64, U::Float64, tol::Float64) = relative_gap(L, U)  > tol || ~(L > -Inf)

#=
"""
    gershgorin_λmin

Computes a lower bound on the smallest eigenvalue of `x` by means of Gershgorin's
Circle Theorem. See the link provided for information
https://mathworld.wolfram.com/GershgorinCircleTheorem.html.
"""
function gershgorin_λmin(x::AbstractMatrix{T}) where T
    xsum = -sum(i -> abs(@inbounds x[i]), diagind(x))
    xmin = typemax(T)
    @inbounds for k in diagind(x)
        xk = x[k]
        xadd = abs(xk) + xk
        xsum += xadd
        if xsum < xmin
            xmin = xsum
        end
        xsum -= xadd
    end
    return xmin
end
=#