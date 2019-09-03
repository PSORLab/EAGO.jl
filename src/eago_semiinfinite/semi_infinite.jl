"""SIP_opts an option type for the SIP routine. The fields it contains are as follows:
- tol: the absolute tolerance for convergence
- kmax: the maximum number of iterations
- eps_g0: initial eps_g restriction parameter, must be greater than zero
- r0: parameter used to reduce eps_g as each iteration (e.g. eps_g = eps_g/r0)
- verbosity: Sets verbosity of the solution routine. Either "None", "Normal", and "Full".
- return_hist: Sets the return type. True returns LBD/UBD and times at each iteration.
- LLP_Solve:
- LBP_Solve:
- UBP_Solve:
- LLP_Opt: LLP solver options.
- LBP_Opt: LBP solver options.
- UBP_Opt: UBP solver options.
- gSIPExp: Expression for semi-infinite inequality constraint
- hSIPExp: Expression for semi-infinite equality constraints
"""
mutable struct SIP_Options
  lower_disc_set::Vector{Vector{Float64}}
  upper_disc_set::Vector{Vector{Float64}}
  tolerance::Float64
  max_iterations::Int
  initial_eps_g::Float64
  initial_r::Float64
  return_hist::Bool
  header_interval::Int
  print_interval::Int
  verbosity::String
  inner_tolerance::Float64
  model
  initialize_bnd_prob
  initialize_extras
end
SIP_Options() = SIP_Options(Vector{Float64}[], Vector{Float64}[], 1E-3, 100000, 1.0,
                            2.0, false, 20, 1, "Normal", 1.0E-8, Model(), nothing, nothing)

"""
    SIP_Result
--------------------------------------------------------------------------------
Description:
Composite type for storing the resulting form SIP solution routine.
--------------------------------------------------------------------------------
Fields:
k             Int64 - Number of iterations run
UBD           Float64 - Upper bound of SIP solution (optimal value)
LBD           Float64 - Lower bound of SIP solution
feas          Bool - Feasibility of SIP
LBP_time      Float64 - Time spent solving the lower bounding problem (sec)
LLP_time      Float64 - Time spent solving the lower level problem (sec)
UBP_time      Float64 - Time spent solving the upper bounding problem (sec)
xbar          Array{Float64} - Solution point
--------------------------------------------------------------------------------
"""
mutable struct SIP_Result
    iteration_number::Int
    upper_bound::Float64
    lower_bound::Float64
    feasibility::Bool
    lower_bounding_time::Float64
    lower_level_time::Float64
    upper_bounding_time::Float64
    x_bar::Vector{Float64}
    p_bar::Vector{Float64}
end
SIP_Result() = SIP_Result(1, Inf, -Inf, true, 0.0, 0.0, 0.0, Float64[], Float64[])

"""
    SIP_Problem_Storage

Storage to pass problem information and solution routine options.
"""
mutable struct SIP_Problem_Storage
  f::Function
  gSIP::Function
  x_l::Vector{Float64}
  x_u::Vector{Float64}
  p_l::Vector{Float64}
  p_u::Vector{Float64}

  np::Int
  nx::Int
  opts::SIP_Options

  h
  hj
  y_l::Vector{Float64}
  y_u::Vector{Float64}
  ny::Int
  upper::Symbol
  sense::Symbol

  xp_l::Vector{Float64}
  xp_u::Vector{Float64}
  n_xp::Int
end


"""
    print_int!(sip_options,k_int,lbd,ubd,eps,r)
--------------------------------------------------------------------------------
Description:
Prints current iteration statistics.
--------------------------------------------------------------------------------
Inputs:
sip_options   SIP_opts - Options object
k_int    Int64   - Iteration Number
lbd      Float64 - Lower Bound
ubd      Float64 - Upper Bound
eps      Float64 - Restriction size
r        Float64 - Restriction Adjustment
--------------------------------------------------------------------------------
Returns:
No returned value but prints the iteration, lower bound, upper bound,
restriction value, r value, absolute ratio
and relative ratio if the verbosity is set to "Normal".
--------------------------------------------------------------------------------
"""
@inline function print_int!(sip_options::SIP_Options, k_int::Int, lbd::Float64, ubd::Float64, eps::Float64, r::Float64)
  if (sip_options.verbosity == "Normal")
    # prints header line every hdr_intv times
    if (mod(k_int,sip_options.header_interval)==0||k_int==1)
      println("Iteration   LBD    UBD     eps      r     Absolute_Gap    Absolute_Ratio")
    end
    # prints iteration summary every prnt_intv times
    if (mod(k_int,sip_options.print_interval)==0)
      ptr_arr_temp = [k_int lbd ubd eps r (ubd-lbd) (lbd/ubd)]
      ptr_arr1 = join([@sprintf("%6u",x) for x in ptr_arr_temp[1]], ",   ")
      ptr_arr2 = join([@sprintf("%3.7f",x) for x in ptr_arr_temp[2:5]], ",     ")
      ptr_arr3 = join([@sprintf("%6u",x) for x in ptr_arr_temp[6:7]], ",")
      println(string(ptr_arr1,",      ",ptr_arr2,",      ",ptr_arr3))
    end
  end
end

function print_llp1!(sip_options::SIP_Options, INNg1::Float64, pbar::Vector{Float64}, feas::Bool)
  if (sip_options.verbosity == "Full" || sip_options.verbosity == "Normal")
    println("solved INN #1: ",INNg1," ",pbar," ",feas)
  end
end

function print_llp2!(sip_options::SIP_Options, INNg2::Float64, pbar::Vector{Float64}, feas::Bool)
  if (sip_options.verbosity == "Full" || sip_options.verbosity == "Normal")
    println("solved INN #2: ",INNg2," ",pbar," ",feas)
  end
end

function print_lbp!(sip_options::SIP_Options, LBDg::Float64, xbar::Vector{Float64}, feas::Bool)
  if (sip_options.verbosity == "Full" || sip_options.verbosity == "Normal")
    println("solved LBD: ",LBDg," ",xbar," ",feas)
  end
end

function print_ubp!(sip_options::SIP_Options, UBD_temp::Float64, xbar::Vector{Float64}, feas::Bool)
  if (sip_options.verbosity == "Full" || sip_options.verbosity == "Normal")
    println("solved UBD: ",UBD_temp," ",xbar," ",feas)
  end
end

function check_inputs(sip_options::SIP_Options)
  if (sip_options.initial_r <= 1.0)
    error("initial_r must be greater than 1")
  elseif (sip_options.initial_eps_g <= 0.0)
    error("eps_g must be greater than 0")
  end
  return sip_options.initial_eps_g, sip_options.initial_r
end

function check_convergence(LBDg::Float64, UBDg::Float64, tolerance::Float64)
  if (abs(UBDg - LBDg) < tolerance)
    println("Algorithm Converged")
    return true
  end
  return false
end

function check_lbp_feasible(feas::Bool)
  if (~feas)
    println("Lower Bounding Problem Not Feasible. Algorithm Terminated")
  end
  return feas
end

function post_llp1_update(INNg1::Float64, inner_tolerance::Float64,
                          xbar::Vector{Float64}, pbar::Vector{Float64},
                          LBDg::Float64, UBDg::Float64, lower_disc_set::Vector{Vector{Float64}})
  xstar = fill(NaN,(length(xbar),))
  UBDg_out = UBDg
  return_flag = false
  if (INNg1 + inner_tolerance <= 0.0)
    xstar[:] = xbar[:]
    UBDg_out = LBDg
    return_flag = true
  else
    push!(lower_disc_set, pbar)
  end
  return xstar, UBDg_out, return_flag
end

function post_llp2_update(INNg2::Float64, inner_tolerance::Float64, r::Float64,
                          eps_g::Float64, xbar::Vector{Float64}, pbar::Vector{Float64},
                          UBD_temp::Float64, UBDg::Float64, upper_disc_set::Vector{Vector{Float64}})
  xstar = fill(NaN,(length(xbar),))
  UBDg_out = UBDg
  eps_g_out = eps_g
  if (INNg2 + inner_tolerance < 0.0)
    if (UBD_temp <= UBDg)
      UBDg_out = UBD_temp
      xstar[:] = xbar[:]
    end
    eps_g_out = eps_g/r
  else
    push!(upper_disc_set, pbar)
  end
  return UBDg_out, xstar, eps_g_out
end

function core_sip_routine(lower_level_problem::Function, bounding_problem::Function,
                          set_xpbar::Function, problem_storage::SIP_Problem_Storage)

  # initializes solution
  UBDg = Inf; LBDg = -Inf; k = 0

  xbar, pbar, nx, np = set_xpbar(problem_storage)
  lower_disc_set = problem_storage.opts.lower_disc_set;
  upper_disc_set = problem_storage.opts.upper_disc_set
  INNg1 = Inf; INNg2 = Inf; feas = true
  xstar = fill(NaN,(nx,))

  sip_result = SIP_Result()
  sip_result.x_bar = xbar
  sip_result.p_bar = pbar
  sip_result.lower_bound = -Inf
  sip_result.upper_bound = Inf

  # checks inputs
  eps_g, r = check_inputs(problem_storage.opts)

  # checks for convergence
  for k=1:problem_storage.opts.max_iterations

    # check for termination
    check_convergence(sip_result.lower_bound, sip_result.upper_bound, problem_storage.opts.tolerance) && (break)

    # solve lower bounding problem and check feasibility
    feas = bounding_problem(lower_disc_set, 0.0, sip_result, problem_storage, true)
    print_lbp!(problem_storage.opts, sip_result.lower_bound, sip_result.x_bar, feas)
    sip_result.feasibility = check_lbp_feasible(feas)
    ~sip_result.feasibility && (break)

    # solve inner program  and update lower discretization set
    INNg1, feas = lower_level_problem(sip_result.x_bar,  sip_result, problem_storage)
    print_llp1!(problem_storage.opts, INNg1, sip_result.p_bar, feas)
    xstar, UBDg, return_flag = post_llp1_update(INNg1, problem_storage.opts.inner_tolerance,
                                                 sip_result.x_bar, sip_result.p_bar, sip_result.lower_bound,
                                                 sip_result.upper_bound, lower_disc_set)
    return_flag && (return sip_result.lower_bound, sip_result.upper_bound, sip_result)

    # solve upper bounding problem, if feasible solve lower level problem,
    # and potentially update upper discretization set
    feas = bounding_problem(upper_disc_set, eps_g, sip_result, problem_storage, false)
    print_ubp!(problem_storage.opts, sip_result.upper_bound, sip_result.x_bar, feas)
    if (feas)
      INNg2, feas = lower_level_problem(sip_result.x_bar,  sip_result, problem_storage)
      print_llp2!(problem_storage.opts, INNg2, sip_result.p_bar, feas)
      UBD_temp = copy(sip_result.upper_bound)
      sip_result.upper_bound, xstar, eps_g = post_llp2_update(INNg2, problem_storage.opts.inner_tolerance,
                                             r, eps_g, sip_result.x_bar, sip_result.p_bar, UBD_temp, UBDg,
                                             upper_disc_set)
    else
      eps_g = eps_g/r
    end

    # print iteration information and advance
    print_int!(problem_storage.opts, k, sip_result.lower_bound, sip_result.upper_bound, eps_g, r)
    sip_result.iteration_number = k
  end

  return sip_result
end

include("sip_explicit.jl")
include("sip_implicit.jl")
