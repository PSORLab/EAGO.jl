module ConvexityRules

    const CV_TRAITS = Dict{Tuple{Union{Expr,Symbol},Symbol,Int},Tuple{Symbol,Symbol,Union{Symbol,Float64},Float64,Float64}}()

    cvtrait(M::Union{Expr,Symbol}, f::Symbol, q::Symbol, args...) = CV_TRAITS[M,f,q,length(args)](args...)

    function _get_quoted_symbol(ex::Expr)
        @assert ex.head == :quote
        @assert length(ex.args) == 1 && isa(ex.args[1], Symbol) "Function not a single symbol"
        ex.args[1]
    end

    function _get_quoted_symbol(ex::QuoteNode)
        @assert isa(ex.value, Symbol) "Function not a single symbol"
        ex.value
    end

    macro set_cvtrait(def,c,m,val,lb,ub)
        @assert isa(def, Expr)
        lhs = def
        @assert isa(lhs, Expr) && lhs.head == :call "LHS is not a dot call"
        qualified_f = lhs.args[1]
        @assert isa(qualified_f, Expr) && qualified_f.head == :(.) "Function is not qualified by module"
        M = qualified_f.args[1]
        f = _get_quoted_symbol(qualified_f.args[2])
        args = lhs.args[2:end]
        key = Expr(:tuple, Expr(:quote, M), Expr(:quote, f), length(args))
        rule = Expr(:tuple, c, m, val, lb, ub)
        return esc(quote
            $ConvexityRules.CV_TRAITS[$key] = $rule
            $key
        end)
    end

    #=
    @set_cvtrait Base.abs2(x)
    @set_cvtrait Base.inv(x)
    =#

    #                                         #Convexity       #Mono     #Switch    #Domain
    @set_cvtrait Base.deg2rad(x)               :Affine      :Increasing    -Inf    -Inf  Inf
    @set_cvtrait Base.rad2deg(x)               :Affine      :Increasing    -Inf    -Inf  Inf
    @set_cvtrait Base.:+(x)                    :Affine      :Increasing    -Inf    -Inf  Inf
    @set_cvtrait Base.transpose(x)             :Affine      :Increasing    -Inf    -Inf  Inf

    @set_cvtrait Base.:-(x)                    :Affine      :Decreasing    -Inf    -Inf  Inf

    @set_cvtrait Base.log(x)                   :Concave     :Increasing    -Inf    0.0   Inf
    @set_cvtrait Base.log10(x)                 :Concave     :Increasing    -Inf    0.0   Inf
    @set_cvtrait Base.log2(x)                  :Concave     :Increasing    -Inf    0.0   Inf
    @set_cvtrait Base.log1p(x)                 :Concave     :Increasing    -Inf   -1.0   Inf
    @set_cvtrait Base.sqrt(x)                  :Concave     :Increasing    -Inf    0.0   Inf
    @set_cvtrait Base.acosh(x)                 :Concave     :Increasing    -Inf    1.0   Inf

    @set_cvtrait Base.exp(x)                   :Convex      :Increasing    -Inf   -Inf   Inf
    @set_cvtrait Base.exp2(x)                  :Convex      :Increasing    -Inf   -Inf   Inf
    @set_cvtrait Base.exp10(x)                 :Convex      :Increasing    -Inf   -Inf   Inf
    @set_cvtrait Base.expm1(x)                 :Convex      :Increasing    -Inf   -Inf   Inf

    @set_cvtrait Base.abs(x)                   :Convex      :DecrToInrc    -Inf   -Inf   Inf
    @set_cvtrait Base.cosh(x)                  :Convex      :DecrToIncr     0.0   -Inf   Inf

    @set_cvtrait Base.atan(x)                  :Convexoconcave :Increasing  0.0  -Inf  Inf
    @set_cvtrait Base.asinh(x)                 :Convexoconcave :Increasing  0.0  -Inf  Inf
    @set_cvtrait Base.tanh(x)                  :Convexoconcave :Increasing  0.0  -Inf  Inf
    @set_cvtrait SpecialFunctions.erf(x)       :Convexoconcave :Increasing  0.0  -Inf  Inf

    @set_cvtrait Base.asech(x)                 :Convexoconcave :Decreasing  0.5   0.0  1.0
    @set_cvtrait Base.acos(x)                  :Convexoconcave :Decreasing  0.0  -1.0  1.0
    @set_cvtrait SpecialFunctions.erfcinv(x)   :Convexoconcave :Decreasing  0.0  -Inf  Inf

    @set_cvtrait Base.tan(x)                   :Concavoconvex  :Increasing  0.0  -pi/2  pi/2
    @set_cvtrait Base.tand(x)                  :Concavoconvex  :Increasing  0.0  -90   90
    @set_cvtrait Base.cbrt(x)                  :Concavoconvex  :Increasing  0.0  -Inf  Inf
    @set_cvtrait Base.sinh(x)                  :Concavoconvex  :Increasing  0.0  -Inf  Inf
    @set_cvtrait Base.asin(x)                  :Concavoconvex  :Increasing  0.0  -1.0  1.0
    @set_cvtrait Base.atanh(x)                 :Concavoconvex  :Increasing  0.0  -1.0  1.0
    @set_cvtrait SpecialFunctions.erfinv(x)    :Concavoconvex  :Increasing  0.0  -1.0  1.0

    @set_cvtrait SpecialFunctions.erfc(x)      :Concavoconvex  :Decreasing  0.0  -Inf  Inf

    @set_cvtrait Base.sech(x)                  :Convexoconcave :IncrToDecr  -Inf -Inf  Inf
end
