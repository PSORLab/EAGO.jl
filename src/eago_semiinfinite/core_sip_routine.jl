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
