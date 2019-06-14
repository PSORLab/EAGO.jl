# defines primitives for CALLUNIVAR operators
for i in (abs,sin, cos, tan, sec, csc, cot, asin, acos, atan, asec, acsc,
          acot, sinh, cosh, tanh, asinh, acosh, atanh, sech, asech, csch,
          acsch, coth, acoth, sqrt, log, log2, log10, log1p, exp, exp2,expm1
          , +, -, inv)
    id = JuMP._Derivatives.univariate_operator_to_id[Symbol(i)]
    @eval function Cassette.overdub(ctx::TraceCtx, ::typeof($i), x::SetTrace)
                node = Tracer.NodeData(JuMP.CALLUNIVAR,$id,[val(x)])
                add_set_node!(ctx.metadata,node)
                return SetTrace(ctx.metadata.set_trace_count)
          end
    @eval Cassette.overdub(ctx::TraceCtx, ::typeof($i), x::Real) = ($i)(x)
end

# defines primitives for bivariate CALL operators (NEED TO ADD ^)
for i in (*,+,-,/,min,max)
    id = JuMP._Derivatives.operator_to_id[Symbol(i)]
    @eval function Cassette.overdub(ctx::TraceCtx, ::typeof($i), x::SetTrace, y::SetTrace)
               node = Tracer.NodeData(JuMP.CALL,$id,[val(x),val(y)])
               add_set_node!(ctx.metadata,node)
               SetTrace(ctx.metadata.set_trace_count)
          end
    for j in (Int16,Int32,Int64,Float16,Float32,Float64)
        @eval function Cassette.overdub(ctx::TraceCtx, ::typeof($i), x::SetTrace, y::($j))
                    c_val = add_constant(ctx.metadata,y)
                    node = Tracer.NodeData(JuMP.CALL,$id,[val(x),c_val])
                    add_set_node!(ctx.metadata,node)
                    SetTrace(ctx.metadata.set_trace_count)
              end
        @eval function Cassette.overdub(ctx::TraceCtx, ::typeof($i), x::($j), y::SetTrace)
                    c_val = add_constant(ctx.metadata,x)
                    node = Tracer.NodeData(JuMP.CALL,$id,[c_val,val(y)])
                    add_set_node!(ctx.metadata,node)
                    SetTrace(ctx.metadata.set_trace_count)
              end
    end
    @eval Cassette.overdub(ctx::TraceCtx, ::typeof($i), x::Real, y::Real) = ($i)(x,y)
end

# defines primitives for bivariate COMPARISON operators
for i in (>,<,==,>=,<=)
    id = JuMP._Derivatives.comparison_operator_to_id[Symbol(i)]
    @eval function Cassette.overdub(ctx::TraceCtx, ::typeof($i), x::SetTrace, y::SetTrace)
                node = Tracer.NodeData(JuMP.COMPARISON,$id,[val(x),val(y)])
                add_set_node!(ctx.metadata,node)
                SetTrace(ctx.metadata.set_trace_count)
            end
    for j in (Int16,Int32,Int64,Float16,Float32,Float64)
        @eval function Cassette.overdub(ctx::TraceCtx, ::typeof($i), x::SetTrace, y::($j))
                    c_val = add_constant(ctx.metadata,y)
                    node = Tracer.NodeData(JuMP.COMPARISON,$id,[val(x),c_val])
                    add_set_node!(ctx.metadata,node)
                    SetTrace(ctx.metadata.set_trace_count)
              end
        @eval function Cassette.overdub(ctx::TraceCtx, ::typeof($i), x::($j), y::SetTrace)
                    c_val = add_constant(ctx.metadata,x)
                    node = Tracer.NodeData(JuMP.COMPARISON,$id,[c_val,val(y)])
                    add_set_node!(ctx.metadata,node)
                    SetTrace(ctx.metadata.set_trace_count)
              end
    end
    @eval Cassette.overdub(ctx::TraceCtx, ::typeof($i), x::Real, y::Real) = ($i)(x,y)
end


# defines primitives overdubbing _if, and _while
#=
function Cassette.overdub(ctx::TraceCtx, ::typeof(_while), x, y, z)
end
function Cassette.overdub(ctx::TraceCtx, ::typeof(_if), x, y)
end
function Cassette.overdub(ctx::TraceCtx, ::typeof(_if), x, y, z)
end
=#
