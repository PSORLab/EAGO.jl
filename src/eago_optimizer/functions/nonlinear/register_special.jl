# Copyright (c) 2018: Matthew Wilhelm, Robert Gottlieb, Dimitri Alston,
# Matthew Stuber, and the University of Connecticut (UConn).
# This code is licensed under the MIT license (see LICENSE.md for full details).
################################################################################
# EAGO
# A development environment for robust and global optimization.
# https://github.com/PSORLab/EAGO.jl
################################################################################
# src/eago_optimizer/functions/nonlinear/register_special.jl
# Defines functions to register nonstandard nonlinear terms in JuMP models.
################################################################################

# Operators with 1 input
for operator in (:relu, :leaky_relu, :maxsig, :maxtanh, :softplus, :pentanh, :sigmoid, :bisigmoid,
                 :softsign, :gelu, :swish, :xabsx, :logcosh, :xlogx, :erf, :erfinv, :erfc)
    func = quote
        model = owner_model(x)
        op = get(model.obj_dict, Symbol("EAGO_", $(String(operator))), nothing)
        if op === nothing
            op = add_nonlinear_operator(
                model,
                1,
                getfield(McCormick, Symbol($(String(operator)))),
                getfield(McCormick, Symbol($(String(operator)), "_deriv")),
                getfield(McCormick, Symbol($(String(operator)), "_deriv2"));
                name = Symbol("EAGO_", $(String(operator))),
            )
            model[Symbol("EAGO_", $(String(operator)))] = op
        end
        return op(x)
    end

    @eval @inline ($operator)(x::VariableRef) = $func
end

# Operator for `McCormick.erfcinv`
function erfcinv(x::VariableRef)
    model = owner_model(x)
    op = get(model.obj_dict, :EAGO_erfcinv, nothing)
    if op === nothing
        op = add_nonlinear_operator(
            model,
            1,
            McCormick.erfcinv,
            x -> -McCormick.erfinv_deriv(1.0 - x),
            x -> McCormick.erfinv_deriv2(1.0 - x);
            name = :EAGO_erfcinv,
        )
        model[:EAGO_erfcinv] = op
    end
    return op(x)
end

# Operators with 2 inputs
for operator in (:xexpax, :arh, :param_relu, :elu)
    func = quote
        model = owner_model(x)
        op = get(model.obj_dict, Symbol("EAGO_", $operator), nothing)
        if op === nothing
            op = add_nonlinear_operator(
                model,
                2,
                getfield(McCormick, Symbol($operator)),
                getfield(McCormick, Symbol($operator, "_grad"));
                name = Symbol("EAGO_", $operator),
            )
            model[Symbol("EAGO_", $operator)] = op
        end
        return op(x, p)
    end

    @eval @inline ($operator)(x::VariableRef, p::VariableRef) = $func
end

# Operator for `McCormick.selu`
function selu(x::VariableRef, α::VariableRef, λ::VariableRef)
    model = owner_model(x)
    op = get(model.obj_dict, :EAGO_selu, nothing)
    if op === nothing
        op = add_nonlinear_operator(
            model,
            3,
            McCormick.selu,
            McCormick.selu_grad;
            name = :EAGO_selu,
        )
        model[:EAGO_selu] = op
    end
    return op(x, α, λ)
end

# Operators for `McCormick.lower_bnd` and `McCormick.upper_bnd`
for operator in (:lower_bnd, :upper_bnd)
    func = quote
        model = owner_model(x)
        op = get(model.obj_dict, Symbol("EAGO_", $operator), nothing)
        if op === nothing
            op = add_nonlinear_operator(
                model,
                2,
                getfield(McCormick, Symbol($operator)),
                getfield(McCormick, Symbol("d_", $operator, "_grad"));
                name = Symbol("EAGO_", $operator),
            )
            model[Symbol("EAGO_", $operator)] = op
        end
        return op(x, b)
    end

    @eval @inline ($operator)(x::VariableRef, b::VariableRef) = $func
end

# Operator for `McCormick.bnd`
function bnd(x::VariableRef, lb::VariableRef, ub::VariableRef)
    model = owner_model(x)
    op = get(model.obj_dict, :EAGO_bnd, nothing)
    if op === nothing
        op = add_nonlinear_operator(
            model,
            3,
            McCormick.bnd,
            McCormick.d_bnd_grad;
            name = :EAGO_bnd,
        )
        model[:EAGO_bnd] = op
    end
    return op(x, lb, ub)
end

# Operators for `McCormick.positive` and `McCormick.negative`
for operator in (:positive, :negative)
    func = quote
        model = owner_model(x)
        op = get(model.obj_dict, Symbol("EAGO_", $operator), nothing)
        if op === nothing
            op = add_nonlinear_operator(
                model,
                1,
                getfield(McCormick, Symbol($operator)),
                x -> 1.0,
                x -> 0.0;
                name = Symbol("EAGO_", $operator),
            )
        model[Symbol("EAGO_", $operator)] = op
        end
        return op(x)
    end

    @eval @inline ($operator)(x::VariableRef) = $func
end
