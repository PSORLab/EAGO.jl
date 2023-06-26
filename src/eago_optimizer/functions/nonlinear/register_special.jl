# Copyright (c) 2018: Matthew Wilhelm, Robert Gottlieb, Dimitri Alston,
# Matthew Stuber, and the University of Connecticut (UConn).
# This code is licensed under the MIT license (see LICENSE.md for full details).
################################################################################
# EAGO
# A development environment for robust and global optimization.
# https://github.com/PSORLab/EAGO.jl
################################################################################
# src/eago_optimizer/functions/nonlinear/register_special.jl
# Defines the function used to register nonstandard nonlinear terms.
################################################################################

"""
register_eago_operators!

Registers all nonstandard nonlinear terms available in EAGO in a JuMP. Uses of
these is generally preferable in EAGO as the relaxations EAGO will generate
will usually be tighter (speeding up convergence time). Note that this will
work can be used by other nonlinear solvers (Ipopt for instance).
"""
function register_eago_operators!(m::JuMP.Model)

    # Register activation functions w/o parameters
    JuMP.register(m, :relu, 1, relu, McCormick.relu_deriv, McCormick.relu_deriv2)
    JuMP.register(m, :leaky_relu, 1, leaky_relu, McCormick.leaky_relu_deriv, McCormick.leaky_relu_deriv2)
    JuMP.register(m, :maxsig, 1, maxsig, McCormick.maxsig_deriv, McCormick.maxsig_deriv2)
    JuMP.register(m, :maxtanh, 1, maxtanh, McCormick.maxtanh_deriv, McCormick.maxtanh_deriv2)
    JuMP.register(m, :softplus, 1, softplus, McCormick.softplus_deriv, McCormick.softplus_deriv2)
    JuMP.register(m, :pentanh, 1, pentanh, McCormick.pentanh_deriv, McCormick.pentanh_deriv2)
    JuMP.register(m, :sigmoid, 1, sigmoid, McCormick.sigmoid_deriv, McCormick.sigmoid_deriv2)
    JuMP.register(m, :bisigmoid, 1, bisigmoid, McCormick.bisigmoid_deriv, McCormick.bisigmoid_deriv2)
    JuMP.register(m, :softsign, 1, softsign, McCormick.softsign_deriv, McCormick.softsign_deriv2)
    JuMP.register(m, :gelu, 1, gelu, McCormick.gelu_deriv, McCormick.gelu_deriv2)
    JuMP.register(m, :swish, 1, swish, McCormick.swish_deriv, McCormick.swish_deriv2)
    JuMP.register(m, :xabsx, 1, xabsx, McCormick.xabsx_deriv, McCormick.xabsx_deriv2)
    JuMP.register(m, :logcosh, 1, logcosh, McCormick.logcosh_deriv, McCormick.logcosh_deriv2)

    # Register activation functions w/ parameters
    MOINL.register_operator(m.nlp_model, :param_relu, 2, param_relu, McCormick.param_relu_grad)
    MOINL.register_operator(m.nlp_model, :elu, 2, elu, McCormick.elu_grad)
    MOINL.register_operator(m.nlp_model, :selu, 3, selu, McCormick.selu_grad)

    # Register other functions
    JuMP.register(m, :xlogx, 1, xlogx, McCormick.xlogx_deriv, McCormick.xlogx_deriv2)
    JuMP.register(m, :f_erf, 1, x -> erf(x), McCormick.erf_deriv, McCormick.erf_deriv2)
    JuMP.register(m, :f_erfinv, 1, x -> erfinv(x), McCormick.erfinv_deriv, McCormick.erfinv_deriv2)
    JuMP.register(m, :f_erfc, 1, x -> erfc(x), McCormick.erfc_deriv, McCormick.erfc_deriv2)
    JuMP.register(m, :f_erfcinv, 1, x -> erfcinv(x), x -> -McCormick.erfinv_deriv(1.0 - x), x -> McCormick.erfinv_deriv2(1.0 - x))

    MOINL.register_operator(m.nlp_model, :arh, 2, arh, McCormick.arh_grad)
    MOINL.register_operator(m.nlp_model, :xexpax, 2, xexpax, McCormick.xexpax_grad)

    # Register bounding functions
    MOINL.register_operator(m.nlp_model, :lower_bnd, 2, lower_bnd, McCormick.d_lower_bnd_grad)
    MOINL.register_operator(m.nlp_model, :upper_bnd, 2, upper_bnd, McCormick.d_upper_bnd_grad)
    MOINL.register_operator(m.nlp_model, :bnd, 3, bnd, McCormick.d_bnd_grad)
    JuMP.register(m, :positive, 1, positive, x -> 1.0, x -> 0.0)
    JuMP.register(m, :negative, 1, negative, x -> 1.0, x -> 0.0)

    return nothing
end
