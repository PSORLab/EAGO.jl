emptyMCinterval(::Type{T}) where T<:AbstractFloat = MCInterval{T}(Inf, -Inf)
emptyMCinterval(x::MCInterval{T}) where T<:AbstractFloat = emptyMCinterval(T)
#const ∅ = emptyMCinterval(Float64)
isempty(x::MCInterval) = x.lo == Inf && x.hi == -Inf

const ∞ = Inf
entireMCinterval(::Type{T}) where T<:AbstractFloat = MCInterval{T}(-Inf, Inf)
entireMCinterval(x::MCInterval{T}) where T<:AbstractFloat = entireMCinterval(T)

function entireinterval(T::Type,V::Type)
    if (V<:MCInterval{T})
        return entireMCinterval(T)
    else
        return IntervalArithmetic.entireinterval(T)
    end
end

isentire(x::MCInterval) = x.lo == -Inf && x.hi == Inf
isunbounded(x::MCInterval) = x.lo == -Inf || x.hi == Inf

isfinite(x::MCInterval) = isfinite(x.lo) && isfinite(x.hi)
isthin(x::MCInterval) = x.lo == x.hi
Base.iszero(x::MCInterval) = iszero(x.lo) && iszero(x.hi)

const float_interval_pi64 = MCInterval(Float64(pi))  # does not change
const float_interval_pi32 = MCInterval(Float32(pi))
const float_interval_pi16 = MCInterval(Float16(pi))
pi_MCinterval(::Type{Float64}) = float_interval_pi64
pi_MCinterval(::Type{Float32}) = float_interval_pi32
pi_MCinterval(::Type{Float16}) = float_interval_pi16
