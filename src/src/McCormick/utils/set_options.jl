"""
    default_options()

Sets McCormickParameters to default.
* `env_max_int = 100`: number of iterations used for envelope calculations
* `env_tol = 1E-10`: tolerance used for envelope calculations
* `mu = 0`: nonsmooth McCormick relaxations used
* `valid_check = true`: relaxtions and intervals are checked for validity
* `subgrad_refine = false`: don't use interval refinement by subgradient propagation
* `multivar_refine = false`: don't use multivariant relaxations
* `mv_tol = 1E-15`: tolerance for multivariant relaxations.
* `outer_rnding = false`: outer rounding of interval and relaxation disabled
* `outer_param = 0.0`: amount of outer rounding
"""
function default_options()
  MC_param = McCormickParamters()
  println("Options set to default values")
end

"""
    set_iterations(val::Integer)

Set iterations for used in envelope calculations.
"""
function set_iterations(val::Integer)
  MC_param.env_max_int = val
  println("Maximum number of iteration for root-finding algorithms used to determing convex/concave envelopes set to $(MC_param.env_max_int)")
end

"""
    set_tolerance(val::Float64)

Set tolerance for used in envelope calculations.
"""
function set_tolerance(val::Float64)
  MC_param.env_tol = val
end

"""
    set_diff_relax(val::Integer)

Set differentiability of relaxations used.
"""
function set_diff_relax(val::Integer)
  diff_relax = val>0
  if (diff_relax>0)
    MC_param.mu = val+1
  else
    MC_param.mu = 0
  end
end

"""
    set_validated(val)

Sets flag for using validity check for McCormick relaxations (cv<=cc).
"""
function set_validated(val)
  MC_param.valid_check = val
end

"""
    set_subgrad_refine(val)

Set flag for using subgradient refinement of interval bounds.
"""
function set_subgrad_refine(val)
  MC_param.subgrad_refine = val
end

"""
    set_multivar_refine(bool,tol)

Set flag for using using multivariant MC relaxations and their tolerance.
"""
function set_multivar_refine(bool,tol)
  MC_param.multivar_refine = bool
  MC_param.mv_tol = tol
end

"""
    set_outer_rnd(bool,tol)

Set flag `bool` for outer round (`cc`/`cv`) and do so by `tol`.
"""
function set_outer_rnd(bool,tol)
  MC_param.outer_rnding = bool
  MC_param.outer_param = tol
end
