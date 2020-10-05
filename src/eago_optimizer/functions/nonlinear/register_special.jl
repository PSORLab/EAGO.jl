

function register_special(r)
    # register activation functions w/o parameters
    register_univariate_operator!(r, :relu, relu, relu_deriv, relu_deriv2)
    register_univariate_operator!(r, :leaky_relu, leaky_relu, leaky_relu_deriv, fprimeprime)
    register_univariate_operator!(r, :maxsig, maxsig, maxsig_deriv, maxsig_deriv2)
    register_univariate_operator!(r, :maxtanh, maxtanh, maxtanh_deriv, maxtanh_deriv2)
    register_univariate_operator!(r, :softplus, softplus, softplus_deriv, softplus_deriv2)
    register_univariate_operator!(r, :pentanh, pentanh, pentanh_deriv, pentanh_deriv2)
    register_univariate_operator!(r, :sigmoid, sigmoid, sigmoid_deriv, fprimeprime)
    register_univariate_operator!(r, :bisigmoid, bisigmoid, bisigmoid_deriv, fprimeprime)
    register_univariate_operator!(r, :softsign, softsign, softsign_deriv, fprimeprime)
    register_univariate_operator!(r, :gelu, gelu, gelu_deriv, fprimeprime)
    register_univariate_operator!(r, :swish1, swish1, swish1_deriv, fprimeprime)

    # register activatio functions w/ parameters
    register_multivariate_operator!(r, :param_relu, f::MOI.AbstractNLPEvaluator)
    register_multivariate_operator!(r, :elu, f::MOI.AbstractNLPEvaluator)
    register_multivariate_operator!(r, :selu, f::MOI.AbstractNLPEvaluator)
end
