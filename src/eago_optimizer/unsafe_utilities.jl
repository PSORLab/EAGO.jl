# Copyright (c) 2018: Matthew Wilhelm & Matthew Stuber.
# This work is licensed under the Creative Commons Attribution-NonCommercial-
# ShareAlike 4.0 International License. To view a copy of this license, visit
# http://creativecommons.org/licenses/by-nc-sa/4.0/ or send a letter to Creative
# Commons, PO Box 1866, Mountain View, CA 94042, USA.
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

Performs `any(f.(x))` in an unsafe manner. Assumes `n == length(x)`. About 4x
faster for small arrays (n < 1000).
"""
function unsafe_any(f, x::Vector{T}, n::Int) where T
    i = 1
    m = n + 1
    while i < m
        v = f(@inbounds x[i])
        if v
            return true
        end
        i += 1
    end
    return false
end

"""
$(TYPEDSIGNATURES)

Performs `any(f.(x, y))` in an unsafe manner. Assumes `n == length(x) == length(y)`.
About 2x faster for small arrays (n < 1000) than `any(x -> (x[1] && x[2]), zip(x, y))`.
"""
function unsafe_any(f, x::Vector{T}, y::Vector{T}, n::Int) where T
    i = 1
    m = n + 1
    while i < m
        v = f(@inbounds x[i], @inbounds y[i])
        if v
            return true
        end
        i += 1
    end
    return false
end


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
