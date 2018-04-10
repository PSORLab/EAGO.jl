# Nonvalidated Interval Library

abstract type AbstractMCInterval{T} <: Real end

struct MCInterval{T<:AbstractFloat} <: AbstractMCInterval{T}
    lo :: T
    hi :: T
end

#MCInterval(a::T, b::T) where T<:AbstractFloat = MCInterval{T}(a, b)
MCInterval(a::T) where T<:AbstractFloat = MCInterval{T}(a, a)
MCInterval(a::Tuple) = MCInterval(a...)
MCInterval(a::T, b::S) where {T<:AbstractFloat, S<:AbstractFloat} = MCInterval(promote(a,b)...)

eltype(x::MCInterval{T}) where T<:AbstractFloat = T
MCInterval(x::MCInterval) = x

# DONE TO HERE!

## FIX CONVERSION ADD CONSTRUCTOR FOR INTERVALS (TO DO!!!!)
#MCInterval{T}(x) where T<:AbstractFloat = MCInterval{T}(convert(T, x))
#MCInterval{T}(x::MCInterval{T}) where T<:AbstractFloat = convert(MCInterval{T}, x)
## FIX CONVERSION ADD CONSTRUCTOR FOR INTERVALS (TO DO!!!!)

include("special.jl")
include("conversion.jl")
#include("precision.jl")
include("set_operations.jl")
include("arithmetic.jl")
include("functions.jl")
include("trigonometric.jl")
include("hyperbolic.jl")
