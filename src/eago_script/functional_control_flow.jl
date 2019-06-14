# THIS IS A DIRECT COPY OF THE FunctionalControlFlow.jl from https://github.com/jekbradbury/FunctionalControlFlow.jl
# We intend to import this module using the package manager if it become available. License below:

#=
MIT License

Copyright (c) 2018 James Bradbury

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
=#

@reexport module FunctionalControlFlow

export @functionalize

using MacroTools: @capture, prewalk, postwalk, isexpr, rmlines, unblock, block

function vars(ex; callables=false, types=false)
    isexpr(ex, Number, QuoteNode) && return Symbol[]
    isexpr(ex, Symbol) && return ex == :nothing ? Symbol[] : [ex]
    ex = ex |> rmlines |> unblock
    @capture(ex, a_.b_) && return vars(a)
    params = Symbol[]
    args = !callables && isexpr(ex, :call) ||
               !types && isexpr(ex, :curly) ? ex.args[2:end] : ex.args
    for x in args
        append!(params, vars(x))
    end
    return params
end
function vars(ex1, ex2; callables=false, types=false)
    params = vars(ex1; callables=callables, types=types)
    append!(params, vars(ex2; callables=callables, types=types))
    unique!(params)
    return params
end

defined(param::Symbol) = :($(Expr(:isdefined, param)) ? $param : nothing)
defined(params::Expr) = Expr(:tuple, (defined(p) for p in params.args)...)

function func(ex, params; mut=true)
    params = :(ret, $(params...))
    ex = block(ex)
    if mut
        ex.args[end] = :(ret = $(ex.args[end]))
        push!(ex.args, :(return $params))
    end
    return Expr(:(->), params, ex)
end
func(ex; mut=true) = func(ex, unique(vars(ex)); mut=mut)

function _while(cond, body, params)
    while cond(params...)
        params = body(params...)
    end
    return params
end

function _if(cond, body, params)
    if cond
        return body(params...)
    end
    return params
end
# TODO consider allowing params1 and params2
function _if(cond, body1, body2, params)
    if cond
        return body1(params...)
    else
        return body2(params...)
    end
end

macro functionalize(ex)
    #TODO break, continue, return
    #continue requires adding conditionals
    #break is continue plus setting an &&-ed while condition to false
    #return from a while is break plus setting a return expr
    postwalk(ex) do x
        isexpr(x, :elseif) && (x.head = :if)
        return if @capture(x, while c_ b_ end |
                              for i_ in v_ b_ end | for i_ = v_ b_ end)
            if isexpr(x, :for)
                c = :(next !== nothing)
                b.args = [:(($i, state) = next), b.args...,
                          :(next = iterate(it, state))]
            end
            params = vars(c, b)
            cf, bf = func(c, params; mut=false), func(b, params)
            ret = quote
                $(bf.args[1]) = _while($cf, $bf, $(defined(bf.args[1])))
                $(bf.args[1].args[1])
            end
            if isexpr(x, :for)
                pushfirst!(ret.args, :(it = $v), :(next = iterate(it)))
            end
            ret
        elseif @capture(x, if c_ b1_ else b2_ end)
            params = vars(b1, b2)
            f1, f2 = func(b1, params), func(b2, params)
            quote
                $(f1.args[1]) = _if($c, $f1, $f2, $(defined(f1.args[1])))
                $(f1.args[1].args[1])
            end
        elseif @capture(x, if c_ b_ end | c_ && b_)
            f = func(b)
            quote
                $(f.args[1]) = _if($c, $f, $(defined(f.args[1])))
                $(f.args[1].args[1])
            end
        elseif @capture(x, a_ && b_)
            f = func(b)
            quote
                c = $a
                $(f.args[1]) = _if(c, $f, $(defined(f.args[1])))
                c ? $(f.args[1].args[1]) : false
            end
        elseif @capture(x, a_ || b_)
            f = func(b)
            quote
                c = !$a
                $(f.args[1]) = _if(c, $f, $(defined(f.args[1])))
                c ? $(f.args[1].args[1]) : true
            end
        elseif isexpr(x, :function)
            x.args[1].args[1] = esc(x.args[1].args[1])
            x
        else
            x
        end
    end
end

end #module
