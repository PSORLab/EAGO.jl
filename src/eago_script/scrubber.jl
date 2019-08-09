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
function generated_scrubbed_method(f::Function, n::Int, inplace = false)
    ci = code_lowered(f)
    @assert(length(ci) == 1)
    function g end
    #if inplace
    #    println("in place gen scub")
        addmethod!(Tuple{typeof(g), Any, Any}, ci[1])
    #else
    #    println("outplace place gen scub")
        addmethod!(Tuple{typeof(g), Any}, ci[1])
    #end
    addmethod!(Tuple{typeof(g),
                     SubArray{Float64,1,Array{Float64,1},Tuple{UnitRange{Int64}},true},
                     SubArray{Float64,1,Array{Float64,1},Tuple{UnitRange{Int64}},true}}, ci[1])
    return g
end

# scrub each udf expression (works for n = 1 )
function scrub(f::Function, n::Int, inplace = false)
    if inplace
        println("ran in place")
        if hasmethod(f, Tuple{Array{Float64,1}})
            println("ran in place 1")
            return false, (y,x) -> overdub(ScrubCtx(), f, y, x)
        else #HERE!!!!
            println("ran in place 2")
            fna = generated_scrubbed_method(f, n, true)
            return true, (y,x) -> overdub(ScrubCtx(), fna, y, x)
        end
    else
        println("ran out place")
        if hasmethod(f, Tuple{Array{Float64,1}})
            println("ran out place 1")
            return false, x -> overdub(ScrubCtx(), f, x)
        else
            println("ran out place 2")
            fna2 = generated_scrubbed_method(f, n)
            return true, x -> overdub(ScrubCtx(), fna2, x)
        end
    end
end

# scrub the model
function scrub!(d::_NLPData)
    # scrub multivariant
    println("start scrub mv")
    user_ops = d.user_operators
    mvop_num = length(user_ops.multivariate_operator_evaluator)
    for i in 1:mvop_num
        evalr = user_ops.multivariate_operator_evaluator[i]
        f_old = evalr.f
        flag, fnew = scrub(evalr.f, evalr.len)
        #ftemp = (out, a) -> gradient!(out, x -> f_old(x...), a)
        #fnew2 = (y, x) -> overdub(ScrubCtx(), ftemp, y, x)
        #flag, fnew2 = scrub(evalr.∇f, evalr.len)
        out = zeros(3)
        a = [1.0; 2.0; 3.0]
        #fnew2(out,a)
        evalr.f = fnew
        #evalr.∇f = fnew2
        user_ops.multivariate_operator_evaluator[i] = evalr
    end
    println("end scrub mv")

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
