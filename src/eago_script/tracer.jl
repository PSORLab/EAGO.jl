# Copyright (c) 2018: Matthew Wilhelm, Robert Gottlieb, Dimitri Alston,
# Matthew Stuber, and the University of Connecticut (UConn).
# This code is licensed under the MIT license (see LICENSE.md for full details).
################################################################################
# EAGO
# A development environment for robust and global optimization.
# https://github.com/PSORLab/EAGO.jl
################################################################################
# src/eago_script/tracer.jl
# Utilities for tracing a user-defined function.
################################################################################

struct NodeInfo
    nodetype::MOINL.NodeType
    index::Int64
    children::Vector{Int64}
end
children(x::NodeInfo) = x.children

# convert method assumes that the children connectivity using in Tracer has
# has been converted to an area hold a single parent value as is used in JuMP
convert(::Type{MOINL.Node}, x::NodeInfo) = MOINL.Node(x.type, x.index, x.children[1])

# val = 1..n corresponds to 1..n variable
struct SetTrace <: Real
    val::Int
end
SetTrace() = SetTrace(0)

val(x::SetTrace) = x.val

struct SetTraceSto
    storage::Vector{SetTrace}
end
getindex(A::SetTraceSto, i::Int64) = getindex(A.storage, i)
iterate(v::SetTraceSto, i=1) = (length(v.storage) < i ? nothing : (v.storage[i], i + 1))

export Tape

# JuMP convention is to store child from function... function call = -1,
# next highest call is has parent 1 and so on... This is a forward tape recording
# architecture, so we'll store children then invert...
mutable struct Tape
    nd::Vector{NodeInfo}
    const_values::Vector{Float64}
    num_valued::Dict{Int,Bool}
    set_trace_count::Int
    const_count::Int
end

Tape() = Tape(NodeInfo[], Float64[], Dict{Int,Bool}(), 0, 0)

function Tape(n::Int)
    node_list = [NodeInfo(MOINL.NODE_VARIABLE,i,[-1]) for i = 1:n]
    num_valued = Dict{Int,Bool}()
    for i = 1:n
        num_valued[i] = false
    end
    Tape(node_list, Float64[], num_valued, n, 0)
end

function add_constant(x::Tape, y)
    x.set_trace_count += 1; x.const_count += 1
    node = NodeInfo(MOINL.NODE_VALUE, x.const_count,[-2]); push!(x.nd, node)
    x.num_valued[x.set_trace_count] = true
    push!(x.const_values, Float64(y))
    x.set_trace_count
end

function add_set_node!(x::Tape, node::NodeInfo)
    push!(x.nd, node)
    x.set_trace_count += 1
    x.num_valued[x.set_trace_count] = false
end


@context TraceCtx

multivariate_operator_to_id = Dict{Symbol,Int}()
for i in 1:length(DEFAULT_MULTIVARIATE_OPERATORS)
    multivariate_operator_to_id[DEFAULT_MULTIVARIATE_OPERATORS[i]] = i
end
univariate_operator_to_id = Dict{Symbol,Int}()
for i in 1:length(DEFAULT_UNIVARIATE_OPERATORS)
    univariate_operator_to_id[DEFAULT_UNIVARIATE_OPERATORS[i]] = i
end
comparison_operators = [:<=, :(==), :>=, :<, :>]
comparison_operator_to_id = Dict{Symbol,Int}()
for i in 1:length(comparison_operators)
    comparison_operator_to_id[comparison_operators[i]] = i
end

# defines primitives for MOINL.NODE_CALL_UNIVARIATE operators
for i in (abs, sin, cos, tan, sec, csc, cot, asin, acos, atan, asec, acsc,
          acot, sinh, cosh, tanh, asinh, acosh, atanh, sech, asech, csch,
          acsch, coth, acoth, sqrt, log, log2, log10, log1p, exp, exp2, expm1,
          +, -, inv)
    id = univariate_operator_to_id[Symbol(i)]
    @eval function Cassette.overdub(ctx::TraceCtx, ::typeof($i), x::SetTrace)
                add_set_node!(ctx.metadata, NodeInfo(MOINL.NODE_CALL_UNIVARIATE, $id, [val(x)]))
                return SetTrace(ctx.metadata.set_trace_count)
          end
    @eval Cassette.overdub(ctx::TraceCtx, ::typeof($i), x::Real) = ($i)(x)
end

# defines primitives for bivariate MOINL.NODE_CALL_MULTIVARIATE operators (NEED TO ADD ^)
for i in (+, -, *, ^, /) # TODO ADD :max, :min
    id = multivariate_operator_to_id[Symbol(i)]
    @eval function Cassette.overdub(ctx::TraceCtx, ::typeof($i), x::SetTrace, y::SetTrace)
               add_set_node!(ctx.metadata, NodeInfo(MOINL.NODE_CALL_MULTIVARIATE, $id, [val(x), val(y)]))
               SetTrace(ctx.metadata.set_trace_count)
          end
    for j in (Int16,Int32,Int64,Float16,Float32,Float64,Irrational)
        @eval function Cassette.overdub(ctx::TraceCtx, ::typeof($i), x::SetTrace, y::($j))
                    add_set_node!(ctx.metadata, NodeInfo(MOINL.NODE_CALL_MULTIVARIATE, $id, [val(x), add_constant(ctx.metadata, y)]))
                    SetTrace(ctx.metadata.set_trace_count)
              end
        @eval function Cassette.overdub(ctx::TraceCtx, ::typeof($i), x::($j), y::SetTrace)
                    add_set_node!(ctx.metadata, NodeInfo(MOINL.NODE_CALL_MULTIVARIATE, $id, [add_constant(ctx.metadata, x),val(y)]))
                    SetTrace(ctx.metadata.set_trace_count)
              end
    end
    @eval Cassette.overdub(ctx::TraceCtx, ::typeof($i), x::Real, y::Real) = ($i)(x,y)
end

# defines primitives for bivariate MOINL.NODE_COMPARISON operators
for i in (>,<,==,>=,<=)
    id = comparison_operator_to_id[Symbol(i)]
    @eval function Cassette.overdub(ctx::TraceCtx, ::typeof($i), x::SetTrace, y::SetTrace)
                add_set_node!(ctx.metadata, NodeInfo(MOINL.NODE_COMPARISON, $id, [val(x), val(y)]))
                SetTrace(ctx.metadata.set_trace_count)
            end
    for j in (Int16,Int32,Int64,Float16,Float32,Float64,Irrational)
        @eval function Cassette.overdub(ctx::TraceCtx, ::typeof($i), x::SetTrace, y::($j))
                    add_set_node!(ctx.metadata, NodeInfo(MOINL.NODE_COMPARISON, $id, [val(x),add_constant(ctx.metadata, y)]))
                    SetTrace(ctx.metadata.set_trace_count)
              end
        @eval function Cassette.overdub(ctx::TraceCtx, ::typeof($i), x::($j), y::SetTrace)
                    add_set_node!(ctx.metadata, NodeInfo(MOINL.NODE_COMPARISON, $id, [add_constant(ctx.metadata, x), val(y)]))
                    SetTrace(ctx.metadata.set_trace_count)
              end
    end
    @eval Cassette.overdub(ctx::TraceCtx, ::typeof($i), x::Real, y::Real) = ($i)(x, y)
end

# define primitives for associative terms
Cassette.overdub(ctx::TraceCtx, ::typeof(*), x, y) = afoldl(x, y)
Cassette.overdub(ctx::TraceCtx, ::typeof(afoldl), x, y, z) = afoldl(x, y, z)
Cassette.overdub(ctx::TraceCtx, ::typeof(afoldl), a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p...) = afoldl(a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p...)

# conversion
Cassette.overdub(ctx::TraceCtx, ::typeof(float), x) = x
Cassette.overdub(ctx::TraceCtx, ::typeof(AbstractFloat), x) = x

# primitive for array access
Cassette.overdub(ctx::TraceCtx, ::typeof(getindex), A::Array, i::Int) = getindex(A,i)
Cassette.overdub(ctx::TraceCtx, ::typeof(getindex), A::SetTraceSto, i::Int) = getindex(A,i)

function Cassette.overdub(ctx::TraceCtx, ::typeof(typeassert), x::Real, type::Type)
    if !isa(x,SetTrace)
        typeassert(x, type)
    end
    return x
end

# prehook for debugging mainly
function Cassette.prehook(ctx::TraceCtx, f::Function, args...)
    #println(f, args)
end

function trace_script(f::Function, n::Int)
    tape = Tape(n)
    if n > 1
        x = SetTraceSto(SetTrace[SetTrace(i) for i = 1:n])
        overdub(TraceCtx(metadata = tape), f, x)
    else
        overdub(TraceCtx(metadata = tape), f, SetTrace(1))
    end
    return tape
end
