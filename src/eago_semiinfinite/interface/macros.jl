# Macro error function
# from InfiniteOpt.jl which was inspired by
# https://github.com/jump-dev/JuMP.jl/blob/709d41b78e56efb4f2c54414266b30932010bd5a/src/macros.jl#L923-L928
function _macro_error(macroname, args, source, str...)
    error("At $(source.file):$(source.line): `@$macroname($(join(args, ", ")))`: ",
          str...)
end

macro decision_variable(model, args...)
    _error(str...) = _macro_error(:decision_variable, (model, args...), __source__, str...)

    extra, kw_args, requestedcontainer = _extract_kw_args(args)
    param_kw_args = filter(kw -> kw.args[1] == :parameter_refs, kw_args)

    # Check for easy case if it is anonymous single variable
    if length(extra) == 0
        code = quote
            @assert isa($model, SIPModel) "Model must be an `SIPModel`."
            JuMP.@variable($model, ($(args...)), variable_type = Infinite, macro_error = $_error)
        end
    else
        x = popfirst!(extra)

        # There are several cases to consider:
        # x                                                     | type of x | x.head
        # ------------------------------------------------------+-----------+------------
        # var                                                   | Symbol    | NA
        # var[1:2]                                              | Expr      | :ref
        # var <= ub or var[1:2] <= ub                           | Expr      | :call
        # lb <= var <= ub or lb <= var[1:2] <= ub               | Expr      | :comparison
        if isexpr(x, :comparison) || isexpr(x, :call)
            inf_expr, params = InfiniteOpt._parse_parameters(_error, Val(x.head),
                                                             x.args)
        else
            inf_expr = x
            params = nothing
        end

        # check for double specification of parameters
        if length(param_kw_args) != 0 && params != nothing
            _error("Cannot specify double specify the infinite parameter references.")
        end

        # easy case where we can don't need to parse extra args
        if length(args) == 1
            code = quote
                @assert isa($model, SIPModel) "Model must be an `SIPModel`."
                JuMP.@variable($model, ($(inf_expr)), variable_type = Infinite, parameter_refs = $params,  macro_error = $_error)
            end
        # here we need to parse the extra args and include them in the call
        else
            rest_args = [args[i] for i = 2:length(args)]
            if isa(params, Nothing)
                code = quote
                    @assert isa($model, SIPModel) "Model must be an `SIPModel`."
                    JuMP.@variable($model, ($(inf_expr)), ($(rest_args...)), variable_type = Infinite, macro_error = $_error)
                end
            else
                code = quote
                    @assert isa($model, SIPModel) "Model must be an `SIPModel`."
                    JuMP.@variable($model, ($(inf_expr)), ($(rest_args...)), variable_type = Infinite, parameter_refs = $params, macro_error = $_error)
                end
            end
        end
    end
    return esc(code)
end
