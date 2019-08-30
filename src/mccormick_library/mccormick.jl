@reexport module McCormick

using DiffRules: diffrule
using StaticArrays: @SVector, SVector, zeros, ones

import Base: +, -, *, /, convert, in, isempty, one, zero, real, eps, max, min,
             abs, inv, exp, exp2, exp10, expm1, log, log2, log10, log1p, acosh, sech,
             csch, coth, acsch, acoth, asech,
             sqrt, sin, cos, tan, min, max, sec, csc, cot, ^, step, sign, intersect,
             promote_rule, asinh, atanh, tanh, atan, asin, cosh, acos,
             sind, cosd, tand, asind, acosd, atand,
             secd, cscd, cotd, asecd, acscd, acotd

import IntervalArithmetic: dist, mid, pow, +, -, *, /, convert, in, isempty,
                           one, zero, real, eps, max, min, abs, exp,
                           expm1, log, log2, log10, log1p, sqrt, ^,
                           sin, cos, tan, min, max, sec, csc, cot, step,sech,
                           csch, coth, acsch, acoth, asech,
                           sign, dist, mid, pow, Interval, interval, sinh, cosh,
                           ∩, IntervalBox, pi_interval, bisect, isdisjoint, length,
                           atan, asin, acos,
                           sind, cosd, tand, asind, acosd, atand,
                           secd, cscd, cotd, asecd, acscd, acotd, half_pi, setrounding

import IntervalContractors: plus_rev, mul_rev, min_rev, max_rev, minus_rev, div_rev, exp_rev,
       exp2_rev, exp10_rev, expm1_rev, log_rev, log2_rev, log10_rev,
       log1p_rev, sin_rev, cos_rev, tan_rev, asin_rev, acos_rev, atan_rev,
       sinh_rev, cosh_rev, tanh_rev, asinh_rev, acosh_rev, atanh_rev,
       abs_rev, sqr_rev, sqrt_rev, power_rev

import Base.MathConstants.golden

# Export forward operators
export MC, cc, cv, Intv, lo, hi,  cc_grad, cv_grad, cnst, +, -, *, /, convert,
       one, zero, dist, real, eps, mid, exp, exp2, exp10, expm1, log, log2,
       log10, log1p, acosh, sqrt, sin, cos, tan, min, max, sec, csc, cot, ^,
       abs, step, sign, pow, in, isempty, intersect, length,
       acos, asin, atan, sinh, cosh, tanh, asinh, atanh, inv, sqr, sech,
       csch, coth, acsch, acoth, asech,
       sind, cosd, tand, asind, acosd, atand,
       sinhd, coshd, tanhd, asinhd, acoshd, atanhd,
       secd, cscd, cotd, asecd, acscd, acotd,
       secdh, cschd, cothd, asechd, acschd, acothd

# Export inplace operators
export plus!, mult!, min!, max!, minus!, div!, exp!, exp2!, exp10!, expm1!,
       log!, log2!, log10!, log1p!, sin!, cos!, tan!, asin!, acos!, atan!,
       sinh!, cosh!, tanh!, asinh!, acosh!, atanh!, abs!, sqr!, sqrt!, pow!

export seed_gradient, IntervalType, set_mc_differentiability!, set_multivar_refine!,
       set_tolerance!, set_iterations!, MC_param

# Export reverse operators
export plus_rev, mul_rev, min_rev, max_rev, minus_rev, div_rev, exp_rev,
       exp2_rev, exp10_rev, expm1_rev, log_rev, log2_rev, log10_rev,
       log1p_rev, sin_rev, cos_rev, tan_rev, asin_rev, acos_rev, atan_rev,
       sinh_rev, cosh_rev, tanh_rev, asinh_rev, acosh_rev, atanh_rev,
       abs_rev, sqr_rev, sqrt_rev, power_rev

# Export utility operators
#=
export grad, zgrad, ∩, mid3, MC_param, mid_grad, seed_g, line_seg, dline_seg,
       outer_rnd, cut, set_valid_check, set_subgrad_refine, set_multivar_refine,
       set_outer_rnd, tighten_subgrad, set_iterations, set_tolerance,
       default_options, value, mincv, maxcc, promote_rule
=#
export mc_opts, gen_expansion_params, gen_expansion_params!, implicit_relax_h,
       implicit_relax_h!, implicit_relax_f, implicit_relax_fg


function __init__()
      setrounding(Interval, :accurate)
end

include("convexity_rules.jl")

include("./mccormick_utilities/constants.jl")
include("mccormick_utilities/inner_utilities.jl")

include("operator_library/type.jl")

include("mccormick_utilities/api_utilities.jl")
include("mccormick_utilities/root_finding.jl")

include("operator_library/forward_operators/forward.jl")
include("operator_library/reverse_operators/reverse.jl")

include("implicit_routines/implicit.jl")

end
