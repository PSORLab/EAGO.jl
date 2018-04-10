# Include files for McCormick operators and export functions

########### exports functions
export SMCg, step, abs, max, min, cos, sin, tan, acos, asin, atan, grad, one, zero
export sqr, pow, inv, sqrt, exp, log, *, +, -, /, ^, cc, cv, lo, hi, convert, dist, real,zgrad
export sinh, cosh, tanh, asinh, acosh, atanh, ∩, mid3, value, mincv, maxcc, promote_rule
export tighten_subgrad, set_iterations, set_tolerance, set_diff_relax, default_options
export set_valid_check, set_subgrad_refine, set_multivar_refine, set_outer_rnd
export MC_param, mid_grad, seed_g, line_seg, dline_seg, outer_rnd, Intv, cut

include("operators/SMCg_Type.jl")
include("utils/utils.jl")
include("utils/root_finding.jl")
include("utils/set_options.jl")

include("operators/SMCg_Power.jl")
include("operators/SMCg_Multiplication.jl")
include("operators/SMCg_Arithmetic.jl")
include("operators/SMCg_ConvexConcave.jl")
include("operators/SMCg_Trignometric.jl")
include("operators/SMCg_Hyperbolic.jl")
include("operators/SMCg_Extrema.jl")
include("operators/SMCg_Conversions.jl")
include("operators/SMCg_Other.jl")

export mc_opts, SetOptions!, MC_KrawczykCW, MC_NewtonGS, GenExpansionParams,
       MC_impRelax, impRelax_f, impRelax_fg, set_default!, InGenExpansionParams,
       MC_NimpRelax, NimpRelax_f, NimpRelax_fg, IndGenExpansionParams,
       MC_NdimpRelax, NdimpRelax_f, NdimpRelax_fg

include("implicit/Options.jl")
include("implicit/Utility.jl")
include("implicit/Contractor.jl")
include("implicit/Affine_Code.jl")
include("implicit/Gen_Param.jl")
include("implicit/Relax_H.jl")
include("implicit/Relax_FG.jl")

mid(x::SMCg) = SMCg(mid(x.Intv),x.IntvBox,x.xref)
