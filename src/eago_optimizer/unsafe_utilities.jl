# Copyright (c) 2018: Matthew Wilhelm & Matthew Stuber.
# This code is licensed under MIT license (see LICENSE.md for full details)
#############################################################################
# EAGO
# A development environment for robust and global optimization
# See https://github.com/PSORLab/EAGO.jl
#############################################################################
# src/eago_optimizer/unsafe_utilities.jl
# Inbounds and non-allocating versions of simple Julia functions.
#############################################################################

"""
$(TYPEDSIGNATURES)

Performs `map!(f, y, x)` in an unsafe manner if y[i] is true, else no-op.
Assumes `n == length(x) == length(y)`. About 2x faster for small arrays (n < 1000).
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
