const IntervalType = Interval{Float64}

const one_intv = one(Interval{Float64})
const half_intv = Interval{Float64}(0.5)
const two_intv = Interval{Float64}(2.0)
const log2_intv = log(Interval{Float64}(2.0))
const log10_intv = log(Interval{Float64}(10.0))

fast_inv(x::Interval{Float64}) = one_intv/x
fast_exp2(x::Interval{Float64}) = exp(log2_intv*x)
fast_exp10(x::Interval{Float64}) = exp(log10_intv*x)

function fast_tanh(x::Interval{Float64})
    t1 = exp(-two_intv*x)
    t2 = exp(two_intv*x)
    pos = (one_intv-t1)/(one_intv+t1)
    neg = (t2-one_intv)/(t2+one_intv)
    return pos ∩ neg
end

function fast_asinh(x::Interval{Float64})
    isempty(x) && return x
    Interval(log(interval(x.lo)+sqrt(one_intv+pow(interval(x.lo),2))).lo,
             log(interval(x.hi)+sqrt(one_intv+pow(interval(x.hi),2))).hi)
end

# ADD Handling for wrong domain, empty domain
function fast_acosh(x::Interval{Float64})
    isempty(x) && return x
    (x.hi <= 1.0) && return Interval(-∞, log(interval(x.hi)+sqrt(interval(x.hi)+one_intv)*sqrt(interval(x.hi)-one_intv)).hi)
    Interval(log(interval(x.lo)+sqrt(interval(x.lo)+one_intv)*sqrt(interval(x.lo)-one_intv)).lo,
             log(interval(x.hi)+sqrt(interval(x.hi)+one_intv)*sqrt(interval(x.hi)-one_intv)).hi)
end

# ADD Handling for wrong domain, empty domain
#=
function fast_atanh(x::Interval{Float64})
    isempty(x) && return x
    (x.hi >= 1.0) && return Interval(half_intv*(log(one_intv+interval(x.lo))-log(one_intv-interval(x.lo))),∞)
    (x.lo <= -1.0)  && return Interval(-∞, half_intv*(log(one_intv+interval(x.hi))-log(one_intv-interval(x.hi))))
    Interval(half_intv*(log(one_intv+interval(x.lo))-log(one_intv-interval(x.lo))),
             half_intv*(log(one_intv+interval(x.hi))-log(one_intv-interval(x.hi))))
end
=#
function fast_pow(x::Interval{Float64},y::Int)
    if y == 0
        return one(Interval{Float64})
    elseif y == 1
        return x
    elseif y > 1
        return pow(x,y)
    else
        return inv(pow(x,-y))
    end
end


# Redefines interval extensions to faster non-allocating versions
@warn("BEGIN REDEFINING INTERVAL OPERATORS TO NON-IEEE-1788-2015 BUT FASTER VERSIONS")
^(x::Interval{Float64},y::Int) = fast_pow(x,y)
^(x::Interval{Float64},y::Float64) = pow(x,y)
# tan(x::Interval{Float64}) = fast_tan(x)
inv(x::Interval{Float64}) = fast_inv(x)
exp2(x::Interval{Float64}) = fast_exp2(x)
exp10(x::Interval{Float64}) = fast_exp10(x)
tanh(x::Interval{Float64}) = fast_tanh(x)
asinh(x::Interval{Float64}) = fast_asinh(x)
acosh(x::Interval{Float64}) = fast_acosh(x)
#atanh(x::Interval{Float64}) = fast_atanh(x)
@warn("END REDEFINING INTERVAL OPERATORS TO NON-IEEE-1788-2015 BUT FASTER VERSIONS")

lo(x::Interval{Float64}) = x.lo
hi(x::Interval{Float64}) = x.hi

function step(x::Interval{Float64})
      isempty(x) && return emptyinterval(x)
      xmin::Float64 = ((x.lo) < 0.0) ? 0.0 : 1.0
      xmax::Float64 = ((x.hi) >= 0.0) ? 1.0 : 0.0
      return Interval{Float64}(xmin,xmax)
end
