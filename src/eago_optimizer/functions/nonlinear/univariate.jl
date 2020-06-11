# Copyright (c) 2018: Matthew Wilhelm & Matthew Stuber.
# This work is licensed under the Creative Commons Attribution-NonCommercial-
# ShareAlike 4.0 International License. To view a copy of this license, visit
# http://creativecommons.org/licenses/by-nc-sa/4.0/ or send a letter to Creative
# Commons, PO Box 1866, Mountain View, CA 94042, USA.
#############################################################################
# EAGO
# A development environment for robust and global optimization
# See https://github.com/PSORLab/EAGO.jl
#############################################################################
# src/eago_optimizer/functions/nonlinear/univariate.jl
# Defines switching functions for forward and reverse univariate functions.
#############################################################################

exprs = Expr[]
for i = 1:length(univariate_operators)
    op = univariate_operators[i]
    ex = :(return $op(x))
    push!(exprs, ex)
end

function binaryswitch(ids, exprs)
    if length(exprs) <= 3
        out = Expr(:if, Expr(:call, :(==), :operator_id, ids[1]), exprs[1])
        if length(exprs) > 1
            push!(out.args, binaryswitch(ids[2:end], exprs[2:end]))
        end
        return out
    else
        mid = length(exprs) >>> 1
        return Expr(:if, Expr(:call, :(<=), :operator_id, ids[mid]),
            binaryswitch(ids[1:mid], exprs[1:mid]),
            binaryswitch(ids[mid+1:end], exprs[mid+1:end]))
    end
end

switchexpr = binaryswitch(1:length(exprs), exprs)
@eval @inline function eval_univariate_set(operator_id,x::T) where T
    $switchexpr
    error("No match for operator_id")
end

#=
null_interval_rev(y::Interval{T}, x::Interval{T}) where T = y, x
const DEFINED_REVERSES = [:plus_rev, :minus_rev, :inv_rev,
                          :mul_rev, :div_rev, power_rev,
                          :sign_rev, max_rev, min_rev,
                          :sqr_rev, sqrt_rev, abs_rev,
                          :exp_rev, exp2_rev, exp10_rev, expm1_rev,
                          :log_rev, log2_rev, log10_rev, log1p_rev,
                          :sin_rev, cos_rev, tan_rev,
                          :asin_rev, acos_rev, atan_rev,
                          :sinh_rev, cosh_rev, tanh_rev,
                          :asinh_rev, acosh_rev, atanh_rev]
                          =#

# code for defining reverse operators
univariate_operators_rev = [:plus_rev, :minus_rev]
for i in 3:length(univariate_operators)
    string = String(univariate_operators[i])*"_rev"
#    if string == "IntervalContractors."
    @show string
    push!(univariate_operators_rev, Symbol(string))
end

univariate_reverse_operator_to_id = Dict{Symbol,Int}()
for i = 1:length(univariate_operators)
    univariate_reverse_operator_to_id[univariate_operators_rev[i]] = i
end

exprs_rev = Expr[]
for i = 1:length(univariate_operators_rev)
    op = univariate_operators_rev[i]
    ex = :(return $op(y,x))
    push!(exprs_rev, ex)
end

switchexpr_rev = binaryswitch(1:length(exprs_rev), exprs_rev)
@eval @inline function eval_univariate_set_reverse(operator_id, y::T, x::T) where T
    $switchexpr_rev
    error("No match for operator_id")
end

function make_uv_gate(symbol_list)
   n = Expr(:||)
   for (i, sym) in enumerate(symbol_list)
      if i !== 1
         n = :(x == $(univariate_operator_to_id[sym]) || $n)
      else
        n = Expr(:||)
        n = :(x == $(univariate_operator_to_id[sym]))
      end
   end
   n
end

single_tp_ops = [:tan, :cot, :asin, :acos, :atan, :sinh, :tanh,
                 :erf, :erfinv, :erfc, :erfcinv, :erfi, :erfcx, :tand, :cotd,
                 :acosd, :atand, :asind, :asinh, :atanh, :asech]
single_tp_gate = make_uv_gate(single_tp_ops)
@eval @inline single_tp(x::Int) = $single_tp_gate

double_tp_ops = [:sin, :cos, :sind, :cosd, :sech, :dawson]
double_tp_gate = make_uv_gate(double_tp_ops)
@eval @inline double_tp(x::Int) = $double_tp_gate

function make_tp_gate_1(symbol_list, x::Symbol, z::Symbol)
    n = Expr(:block)
    for (i, sym) in enumerate(symbol_list)
        a = :(op == $(univariate_operator_to_id[sym]) && (return (McCormick.$sym)(x), 0.0, 0.0))
        push!(n.args, a)
    end
    push!(n.args, :(error("No operator")))

    nd = Expr(:block)
    for (i, sym) in enumerate(symbol_list)
        sym_kernel = Symbol(String(sym)*"_kernel")
        a = :(op == $(univariate_operator_to_id[sym]) && (return (McCormick.$sym_kernel)(x, z.Intv, tp1, tp2)))
        push!(nd.args, a)
    end
    push!(nd.args, :(error("No operator")))

    qn = quote
            if flag
                $n
            else
                $nd
            end
        end
    qn
end

function make_tp_gate_2(symbol_list, x::Symbol, z::Symbol)
    n = Expr(:block)
    for (i, sym) in enumerate(symbol_list)
        a = :(op == $(univariate_operator_to_id[sym]) && (return (McCormick.$sym)(x), 0.0, 0.0, 0.0, 0.0))
        push!(n.args, a)
    end
    push!(n.args, :(error("No operator")))

    nd = Expr(:block)
    for (i, sym) in enumerate(symbol_list)
        sym_kernel = Symbol(String(sym)*"_kernel")
        a = :(op == $(univariate_operator_to_id[sym]) && (return (McCormick.$sym_kernel)(x, z.Intv, tp1, tp2, tp3, tp4)))
        push!(nd.args, a)
    end
    push!(nd.args, :(error("No operator")))

    qn = quote
            if flag
                $n
            else
                $nd
            end
        end
    qn
end

single_tp_expr = make_tp_gate_1(single_tp_ops, :x, :z)
@eval @inline function single_tp_set(op::Int64, x::MC, z::MC, tp1::Float64,
                               tp2::Float64, flag::Bool)
    $single_tp_expr
end
double_tp_expr = make_tp_gate_2(double_tp_ops, :x, :z)
@eval @inline function double_tp_set(op::Int64, x::MC, z::MC, tp1::Float64,
                               tp2::Float64, tp3::Float64, tp4::Float64, flag::Bool)
    $double_tp_expr
end
