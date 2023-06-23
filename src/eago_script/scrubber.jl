# Copyright (c) 2018: Matthew Wilhelm, Robert Gottlieb, Dimitri Alston,
# Matthew Stuber, and the University of Connecticut (UConn).
# This code is licensed under the MIT license (see LICENSE.md for full details).
################################################################################
# EAGO
# A development environment for robust and global optimization.
# https://github.com/PSORLab/EAGO.jl
################################################################################
# src/eago_script/scrubber.jl
# A context used to "scrub" type specific storage, assertions, and other
# language features from a user-defined function.
################################################################################

@context ScrubCtx

# Cassette specific functions for scrubbing udfs of objects that interfere with overloading
Cassette.overdub(ctx::ScrubCtx, ::typeof(typeassert), x::Real, type::Type) = x
function Cassette.overdub(ctx::ScrubCtx, ::typeof(zeros), t, dims...)
    if t <: AbstractFloat
        return zeros(Real, dims...)
    elseif t <: Integer
        return zeros(Real, dims...)
    end
    return zeros(t, dims...)
end
#overdub(ctx::ScrubCtx, ::typeof(zero), type::Type) = zero(type)
function Cassette.overdub(ctx::ScrubCtx, ::typeof(hcat), A...)
    vA = hcat(A...)
    sz = size(vA)
    vR = zeros(Real, sz...)
    vR[:] = vA[:]
    return vR
end
function Cassette.overdub(ctx::ScrubCtx, ::typeof(vcat), A...)
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
    addmethod!(Tuple{typeof(g), Any, Any}, ci[1])
    addmethod!(Tuple{typeof(g), Any}, ci[1])
    addmethod!(Tuple{typeof(g),
                     SubArray{Float64,1,Array{Float64,1},Tuple{UnitRange{Int64}},true},
                     SubArray{Float64,1,Array{Float64,1},Tuple{UnitRange{Int64}},true}}, ci[1])
    return g
end

"""
    scrub(f::Function, n::Int, inplace = false)

Replaces storage objects and circumvents assertions that present a UDF from
being overloaded with subtype Real objects by creating a function which
overdubs `f` in `ScrubCtx`.
"""
function scrub(f::Function, n::Int, inplace = false)
    if inplace
        if hasmethod(f, Tuple{Array{Float64,1}})
            println("scrub arc 1...")
            return false, (y,x) -> overdub(ScrubCtx(), f, y, x)
        else
            println("scrub arc 2...")
            fna = generated_scrubbed_method(f, n, true)
            return true, (y,x) -> overdub(ScrubCtx(), fna, y, x)
        end
    else
        if hasmethod(f, Tuple{Array{Float64,1}})
            return false, x -> overdub(ScrubCtx(), f, x)
        else
            fna2 = generated_scrubbed_method(f, n)
            return true, x -> overdub(ScrubCtx(), fna2, x)
        end
    end
end

"""
    scrub!(d::MOI.NLPBlockData)

Applies scrub to every user-defined function in the a `MOI.NLPBlockData` structure.
"""
function scrub!(d::MOI.NLPBlockData)
    error("Function not updated for newest version of EAGO. Please submit an issue.")
    # scrub multivariant
    user_ops = d.user_operators
    mvop_num = length(user_ops.multivariate_operator_evaluator)
    for i in 1:mvop_num
        evalr = user_ops.multivariate_operator_evaluator[i]
        f_old = evalr.f
        flag, fnew = scrub(evalr.f, evalr.len)
        ftemp = (out, a) -> gradient!(out, x -> f_old(x), a)
        fnew2 = (y, x) -> overdub(ScrubCtx(), ftemp, y, x)
        evalr.f = fnew
        evalr.∇f = fnew2
        user_ops.multivariate_operator_evaluator[i] = evalr
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