# Copyright (c) 2018: Matthew Wilhelm & Matthew Stuber.
# This code is licensed under MIT license (see LICENSE.md for full details)
#############################################################################
# EAGO
# A development environment for robust and global optimization
# See https://github.com/PSORLab/EAGO.jl
#############################################################################
# Defines the a function used to register nonstandard nonlinear terms.
#############################################################################

"""
register_eago_operators!

Registers all nonstandard nonlinear terms available in EAGO in a JuMP. Uses of
these is generally preferable in EAGO as the relaxations EAGO will generate
will usually be tighter (speeding up convergence time). Note that this will
work can be used by other nonlinear solvers (Ipopt for instance).
"""
function register_eago_operators!(m::JuMP.Model)

    # register activation functions w/o parameters
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
    JuMP.register(m, :swish1, 1, swish1, McCormick.swish1_deriv, McCormick.swish1_deriv2)

    # register activatio functions w/ parameters
    d_param_relu = JuMP._UserFunctionEvaluator(x -> param_relu(x...),
                                               (g,x) -> McCormick.param_relu_grad(g,x...), 2)
    d_elu = JuMP._UserFunctionEvaluator(x -> elu(x...), (g,x) -> McCormick.elu_grad(g,x...), 2)
    d_selu = JuMP._UserFunctionEvaluator(x -> selu(x...), (g,x) -> McCormick.selu_grad(g,x...), 3)
    JuMP._Derivatives.register_multivariate_operator!(m.nlp_data.user_operators, :param_relu, d_param_relu)
    JuMP._Derivatives.register_multivariate_operator!(m.nlp_data.user_operators, :elu, d_elu)
    JuMP._Derivatives.register_multivariate_operator!(m.nlp_data.user_operators, :selu, d_selu)

    # register other functions
    JuMP.register(m, :xlogx, 1, xlogx, McCormick.xlogx_deriv, McCormick.xlogx_deriv2)
    JuMP.register(m, :f_erf, 1, x -> erf(x), McCormick.erf_deriv, McCormick.erf_deriv2)
    JuMP.register(m, :f_erfinv, 1, x -> erfinv(x), McCormick.erfinv_deriv, McCormick.erfinv_deriv2)
    JuMP.register(m, :f_erfc, 1, x -> erfc(x), McCormick.erfc_deriv, McCormick.erfc_deriv2)
    JuMP.register(m, :f_erfcinv, 1, x -> erfcinv(x), x -> -McCormick.erfinv_deriv(1.0 - x), x -> McCormick.erfinv_deriv2(1.0 - x))

    d_arh = JuMP._UserFunctionEvaluator(x -> McCormick.arh(x...), (g,x) -> McCormick.arh_grad(g,x...), 2)
    d_xexpax = JuMP._UserFunctionEvaluator(x -> McCormick.xexpax(x...), (g,x) -> McCormick.xexpax_grad(g,x...), 2)
    JuMP._Derivatives.register_multivariate_operator!(m.nlp_data.user_operators, :arh, d_arh)
    JuMP._Derivatives.register_multivariate_operator!(m.nlp_data.user_operators, :xexpax, d_xexpax)

    # register bounding functions
    d_lower_bnd = JuMP._UserFunctionEvaluator(x -> lower_bnd(x...), (g,x) -> McCormick.d_lower_bnd_grad(g,x...), 2)
    d_upper_bnd = JuMP._UserFunctionEvaluator(x -> upper_bnd(x...), (g,x) -> McCormick.d_upper_bnd_grad(g,x...), 2)
    d_bnd = JuMP._UserFunctionEvaluator(x -> bnd(x...), (g,x) -> McCormick.d_bnd_grad(g,x...), 3)
    JuMP._Derivatives.register_multivariate_operator!(m.nlp_data.user_operators, :lower_bnd, d_lower_bnd)
    JuMP._Derivatives.register_multivariate_operator!(m.nlp_data.user_operators, :upper_bnd, d_upper_bnd)
    JuMP._Derivatives.register_multivariate_operator!(m.nlp_data.user_operators, :bnd, d_bnd)
    JuMP.register(m, :positive, 1, positive, x -> 1.0, x -> 0.0)
    JuMP.register(m, :negative, 1, negative, x -> 1.0, x -> 0.0)

    m.nlp_data.largest_user_input_dimension = 3

    return nothing
end
