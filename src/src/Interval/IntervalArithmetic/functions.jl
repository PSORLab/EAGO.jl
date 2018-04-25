# Integer power:
function ^(a::MCInterval{T}, n::Integer) where {T<:AbstractFloat}
    isempty(a) && return a
    n == 0 && return one(a)
    n == 1 && return a
    n == 2 && return sqr(a)
    n < 0 && a == zero(a) && return emptyMCinterval(T)

    if isodd(n) # odd power
        isentire(a) && return a
        if n > 0
            a.lo == 0 && return MCInterval{T}(zero(T), a.hi^n)
            a.hi == 0 && return MCInterval{T}(a.lo^n, zero(T))
            return MCInterval{T}(a.lo^n, a.hi^n)
        else
            if a.lo ≥ 0
                a.lo == 0 && return MCInterval{T}(a.hi^n, infty(T))
                return MCInterval{T}(a.hi^n, a.lo^n)
            elseif a.hi ≤ 0

                a.hi == 0 && return MCInterval{T}(ninfty(T), a.lo^n)
                return MCInterval{T}(a.hi^n, a.lo^n)
            else
                return entireMCinterval(T)
            end
        end

    else # even power
        if n > 0
            if a.lo ≥ 0
                return MCInterval{T}(a.lo^n, a.hi^n)
            elseif a.hi ≤ 0
                return MCInterval{T}(a.hi^n, a.lo^n)
            else
                return MCInterval{T}(mig(a)^n, mag(a)^n)
            end

        else
            if a.lo ≥ 0
                return MCInterval{T}(a.hi^n, a.lo^n)
            elseif a.hi ≤ 0
                return MCInterval{T}(a.lo^n, a.hi^n)
            else
                return MCInterval{T}(mag(a)^n, mig(a)^n)
            end
        end
    end
end

function sqr(a::MCInterval{T}) where {T<:AbstractFloat}
    return a*a
end

function sqrt(a::MCInterval{T}) where {T<:AbstractFloat}
    domain = MCInterval{T}(zero(T), Inf)
    a = a ∩ domain

    isempty(a) && return a

    MCInterval{T}(sqrt(a.lo), sqrt(a.hi))  # `sqrt` is correctly-rounded
end

function pow(x::MCInterval{T}, n::Q) where {T<:AbstractFloat,Q<:Integer}  # fast integer power

    isempty(x) && return x

    if iseven(n) && zero(T) ∈ x
        return hull(zero(x),
                    hull(Base.power_by_squaring(MCInterval(mig(x)), n),
                        Base.power_by_squaring(MCInterval(mag(x)), n)))
    else
        Base.power_by_squaring(MCInterval(x.lo), n)
        Base.power_by_squaring(MCInterval(x.hi), n)
        hull( Base.power_by_squaring(MCInterval(x.lo), n),
                      Base.power_by_squaring(MCInterval(x.hi), n) )
      return hull( Base.power_by_squaring(MCInterval(x.lo), n),
                    Base.power_by_squaring(MCInterval(x.hi), n) )
    end
end

function pow(x::MCInterval, y::Real)  # fast real power, including for y an Interval

    isempty(x) && return x

    return exp(y * log(x))
end




for f in (:exp, :expm1, :exp2, :exp10)
    @eval begin
        function ($f)(a::MCInterval{T}) where T
            isempty(a) && return a
            MCInterval{T}(($f)(a.lo), ($f)(a.hi))
        end
    end
end


for f in (:log, :log2, :log10, :log1p)

    @eval function ($f)(a::MCInterval{T}) where T
            domain = MCInterval{T}(zero(T), Inf)
            a = a ∩ domain

            (isempty(a) || a.hi ≤ zero(T)) && return emptyMCinterval(a)

            MCInterval{T}(($f)(a.lo), ($f)(a.hi))

        end
end
