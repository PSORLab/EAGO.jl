using BenchmarkTools, EAGO.McCormick, IntervalArithmetic

# don't really want to measure interval performance, so we're using the 
# fastest interval rounding mode
setrounding(Interval, :none)

const SUITE = BenchmarkGroup()

S = SUITE["Constructors"] = BenchmarkGroup()
S["MC(Float64)"] = @benchmarkable MC{5,NS}(2.1)
S["MC(Float64,Float64)"] = @benchmarkable MC{5,NS}(0.1,2.0)
S["MC(Interval{Float64})"] = @benchmarkable MC{5,NS}($(Interval{Float64}(0.1,2.0)))
S["MC(Float64, Interval{Float64}, Int)"] = @benchmarkable MC{5,NS}(1.0,$(Interval{Float64}(0.1,2.0)),2)

S = SUITE["OOP NS McCormick"] = BenchmarkGroup()
begin
    for op in (+, -, exp, exp2, exp10, expm1, log, log2, log10, log1p,
               abs, one, zero, real, max, min, inv, cosh, acosh, sqrt,
               isone, isnan)
        S[string(op)*"(x)"] = @benchmarkable $(op)(a) setup = (a = MC{5,NS}(1.0,$(Interval{Float64}(0.1,2.0)),2))
    end
    for op in (asin, sin, asind, sind, acos, acosd, cos, cosd, atan, atand,
               tan, tand, sec, secd, step, sign, sinh, atanh, tanh, )
        S[string(op)*"(x), x > 0"] = @benchmarkable $(op)(a) setup = (a = MC{5,NS}(0.4,$(Interval{Float64}(0.1,0.9)),2))
        S[string(op)*"(x), 0 ∈ x"] = @benchmarkable $(op)(a) setup = (a = MC{5,NS}(0.1,$(Interval{Float64}(-0.4,0.3)),2))
        S[string(op)*"(x), x < 0"] = @benchmarkable $(op)(a) setup = (a = MC{5,NS}(-0.5,$(Interval{Float64}(-0.9,-0.1)),2))
    end
    for op in (min, max, *, -, +, /)
        S[string(op)*"(x, Float64)"] =  @benchmarkable $(op)(x,q) setup = (x = MC{5,NS}(0.4,$(Interval{Float64}(0.1,0.9)),2); q = 1.34534)
        S[string(op)*"(Float64, x)"] =   @benchmarkable $(op)(q,x) setup = (x = MC{5,NS}(0.4,$(Interval{Float64}(0.1,0.9)),2); q = 1.34534)
        S[string(op)*"(x, y)"] =   @benchmarkable $(op)(x,y) setup = (x = MC{5,NS}(0.4,$(Interval{Float64}(0.1,0.9)),2);
                                                                        y = MC{5,NS}(0.5,$(Interval{Float64}(0.3,0.9)),2);)
    end
    S["x^2, x > 0"] = @benchmarkable a^2 setup = (a = MC{5,NS}(0.4,$(Interval{Float64}(0.1,0.9)),2))
    S["x^2, 0 ∈ x"] = @benchmarkable a^2 setup = (a = MC{5,NS}(0.1,$(Interval{Float64}(-0.4,0.3)),2))
    S["x^2, x < 0"] = @benchmarkable a^2 setup = (a = MC{5,NS}(-0.5,$(Interval{Float64}(-0.9,-0.1)),2))
    S["x^3, x > 0"] = @benchmarkable a^3 setup = (a = MC{5,NS}(0.4,$(Interval{Float64}(0.1,0.9)),2))
    S["x^3, 0 ∈ x"] = @benchmarkable a^3 setup = (a = MC{5,NS}(0.1,$(Interval{Float64}(-0.4,0.3)),2))
    S["x^3, x < 0"] = @benchmarkable a^3 setup = (a = MC{5,NS}(-0.5,$(Interval{Float64}(-0.9,-0.1)),2))
    S["x^4, x > 0"] = @benchmarkable a^4 setup = (a = MC{5,NS}(0.4,$(Interval{Float64}(0.1,0.9)),2))
    S["x^4, 0 ∈ x"] = @benchmarkable a^4 setup = (a = MC{5,NS}(0.1,$(Interval{Float64}(-0.4,0.3)),2))
    S["x^4, x < 0"] = @benchmarkable a^4 setup = (a = MC{5,NS}(-0.5,$(Interval{Float64}(-0.9,-0.1)),2))
    S["x^-1, x > 0"] = @benchmarkable a^-1 setup = (a = MC{5,NS}(0.4,$(Interval{Float64}(0.1,0.9)),2))
    S["x^-1, x < 0"] = @benchmarkable a^-1 setup = (a = MC{5,NS}(-0.5,$(Interval{Float64}(-0.9,-0.1)),2))
    S["x^-2, x > 0"] = @benchmarkable a^-2 setup = (a = MC{5,NS}(0.4,$(Interval{Float64}(0.1,0.9)),2))
    S["x^-2, x < 0"] = @benchmarkable a^-2 setup = (a = MC{5,NS}(-0.5,$(Interval{Float64}(-0.9,-0.1)),2))
    S["x^-3, x > 0"] = @benchmarkable a^-3 setup = (a = MC{5,NS}(0.4,$(Interval{Float64}(0.1,0.9)),2))
    S["x^-3, x < 0"] = @benchmarkable a^-3 setup = (a = MC{5,NS}(-0.5,$(Interval{Float64}(-0.9,-0.1)),2))
end
