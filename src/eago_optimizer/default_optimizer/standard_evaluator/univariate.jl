# code for defining forward operators
const univariate_operators = Symbol[:+,:-,:abs]
const univariate_operator_to_id = Dict{Symbol,Int}(:+ => 1, :- => 2, :abs => 3)

for (op, deriv) in Calculus.symbolic_derivatives_1arg()
    push!(univariate_operators, op)
    univariate_operator_to_id[op] = length(univariate_operators)
end

exprs = Expr[]
for i = 1:length(univariate_operators)
    op = univariate_operators[i]
    ex = :(return $op(x))
    push!(exprs, ex)
end

function binaryswitch(ids, exprs)
    #println("exprs: $exprs")
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

# code for defining reverse operators
univariate_operators_rev = [:plus_rev, :minus_rev]
for i in 3:length(univariate_operators)
    string = String(univariate_operators[i])*"_rev"
    push!(univariate_operators_rev, Symbol(string))
end

univariate_reverse_operator_to_id = Dict{Symbol,Int}()
for i in 1:length(univariate_operators)
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
