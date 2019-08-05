@context ScrubCtx

# Cassette specific functions for scrubbing udfs of objects that interfere with overloading
overdub(ctx::ScrubCtx, ::typeof(typeassert), x::Real, type::Type) = x
overdub(ctx::ScrubCtx, ::typeof(zeros), t, dims...) = zeros(Real, dims...)
function overdub(ctx::ScrubCtx, ::typeof(hcat), A...)
    vA = hcat(A...)
    sz = size(vA)
    vR = zeros(Real, sz...)
    vR[:] = vA[:]
    return vR
end
function overdub(ctx::ScrubCtx, ::typeof(vcat), A...)
    vA = vcat(A...)
    sz = size(vA)
    vR = zeros(Real, sz...)
    vR[:] = vA[:]
    return vR
end

# scrub type assertion from definition of udf
function generated_scrubbed_method(f::Function, n::Int)
    ci = code_lowered(f)
    @assert(length(ci) == 1)
    function g end
    addmethod!(Tuple{typeof(g), Any}, ci[1])
    return g
end

# scrub each udf expression (works for n = 1 )
function scrub(f::Function, n::Int)
    if hasmethod(f, Tuple{Array{Float64,1}})
        return false, x -> overdub(ScrubCtx(), f, x)
    else
        fna = generated_scrubbed_method(f, n)
        return true, x -> overdub(ScrubCtx(), fna, x)
    end
end

# scrub the model
function scrub!(d::_NLPData)
    # scrub multivariant
    user_ops = d.user_operators
    mvop_num = length(user_ops.multivariate_operator_evaluator)
    for i in 1:mvop_num
        eval = user_ops.multivariate_operator_evaluator[i]
        flag, eval.f = scrub(eval.f, eval.len)
        if flag
            cfg = GradientConfig(eval.f, zeros(Float64, eval.len))
            ∇f = (out, y) -> gradient!(out, eval.f, y, cfg)
            eval.∇f = ∇f
        end
        flag, eval.∇f = scrub(eval.∇f, eval.len)
        user_ops.multivariate_operator_evaluator[i] = eval
    end

    # scrub univariants
    svop_num = length(user_ops.univariate_operator_f)
    for i in 1:svop_num
        flag, fnew = scrub(user_ops.univariate_operator_f[i], 1)
        user_ops.univariate_operator_f[i] = fnew
        if flag
            f = user_ops.univariate_operator_f[i]
            fprime = x -> derivative(f, x)
            user_ops.univariate_operator_fprime[i] = fprime
            user_ops.univariate_operator_fprimeprime[i] = x -> derivative(fprime, x)
        end
        fnew = scrub(user_ops.univariate_operator_fprime[i], 1)
        user_ops.univariate_operator_fprime[i] = fnew
        fnew = scrub(user_ops.univariate_operator_fprimeprime[i], 1)
        user_ops.univariate_operator_fprimeprime[i] = fnew
    end

    # reassign evaluator
    d.user_operators = user_ops
end
