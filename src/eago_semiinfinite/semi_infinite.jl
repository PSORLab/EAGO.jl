  mutable struct SIPResult
      iteration_number::Int
      upper_bound::Float64
      lower_bound::Float64
      feasibility::Bool
      xsol::Vector{Float64}
      solution_time::Float64
  end
  SIPResult() = SIPResult(1, Inf, -Inf, true, Float64[], 0.0)

  """
      SIP_Problem_Storage

  Storage to pass problem information and solution routine options.
  """
mutable struct SIPProblem
    x_l::Vector{Float64}
    x_u::Vector{Float64}
    p_l::Vector{Float64}
    p_u::Vector{Float64}

    np::Int
    nSIP::Int
    nx::Int
    sense::Symbol

    init_lower_disc::Vector{Vector{Vector{Float64}}}
    init_upper_disc::Vector{Vector{Vector{Float64}}}

    absolute_tolerance::Float64
    iteration_limit::Int
    initial_eps_g::Float64
    initial_r::Float64

    return_hist::Bool
    header_interval::Int
    print_interval::Int
    verbosity::Int
  end
  function SIPProblem(x_l::Vector{Float64}, x_u::Vector{Float64}, p_l::Vector{Float64}, p_u::Vector{Float64}, kwargs)

      np = length(p_l)
      nx = length(x_l)

      absolute_tolerance = haskey(kwargs, :absolute_tolerance) ? kwargs[:return_hist] : 1E-1
      iteration_limit = haskey(kwargs, :iteration_limit) ? kwargs[:return_hist] : 1000
      initial_eps_g = haskey(kwargs, :initial_eps_g) ? kwargs[:return_hist] : 0.9
      initial_r = haskey(kwargs, :initial_r) ? kwargs[:return_hist] : 1.5

      return_hist = haskey(kwargs, :return_hist) ? kwargs[:return_hist] : false
      header_interval = haskey(kwargs, :header_interval) ? kwargs[:header_interval] : 20
      print_interval = haskey(kwargs, :print_interval) ? kwargs[:print_interval] : 1
      verbosity = haskey(kwargs, :verbosity) ? kwargs[:verbosity] : 1

      sense = haskey(kwargs, :sense) ? kwargs[:sense] : :min
      init_lower_disc = haskey(kwargs, :init_lower_disc) ? kwargs[:init_lower_disc] : Vector{Vector{Float64}}[]
      init_upper_disc = haskey(kwargs, :init_upper_disc) ? kwargs[:init_upper_disc] : Vector{Vector{Float64}}[]

      nSIP = count(x -> String(x)[1:4] === "gSIP", keys(kwargs))

      SIPProblem(x_l, x_u, p_l, p_u, np, nSIP, nx, sense, init_lower_disc,
                 init_upper_disc, absolute_tolerance, iteration_limit,
                 initial_eps_g, initial_r, return_hist, header_interval,
                 print_interval, verbosity)
  end

  function print_int!(verbosity::Int, hdr_intv::Int, prt_intv::Int, k_int::Int,
                      lbd::Float64, ubd::Float64, eps::Float64, r::Float64)

    if (verbosity == 1 || verbosity == 2)

      # prints header line every hdr_intv times
      if (mod(k_int, hdr_intv) == 0 || k_int == 1)
        println("| Iteration | Lower Bound | Upper Bound |   eps   |   r   |  Gap  |  Ratio  |")
      end

      # prints iteration summary every prnt_intv times
      if mod(k_int, prt_intv) == 0

        print_str = "| "

        max_len = 9
        temp_str = string(k_int)
        len_str = length(temp_str)
        print_str *= (" "^(max_len - len_str))*temp_str*" | "

        max_len = 11
        temp_str = string(round(lbd, digits = 6))
        len_str = length(temp_str)
        print_str *= (" "^(max_len - len_str))*temp_str*" | "

        max_len = 11
        temp_str = string(round(ubd, digits = 6))
        len_str = length(temp_str)
        print_str *= (" "^(max_len - len_str))*temp_str*" | "

        max_len = 7
        temp_str = string(round(eps, digits = 6))
        len_str = length(temp_str)
        print_str *= (" "^(max_len - len_str))*temp_str*" | "

        max_len = 5
        temp_str = string(round((ubd-lbd), digits = 5))
        len_str = length(temp_str)
        print_str *= (" "^(max_len - len_str))*temp_str*" | "

        max_len = 7
        temp_str = string(round((ubd-lbd)/abs(ubd), digits = 6))
        len_str = length(temp_str)
        print_str *= (" "^(max_len - len_str))*temp_str*" |"

        println(print_str)
      end
    end
    return
  end

  function print_summary!(verb::Int, val::Float64, x::Vector{Float64}, feas::Bool, desc::String)
    (verb == 2 || verb == 1) && println("solved $desc: ", val, " ", x, " ", feas)
    return
  end

  function check_inputs!(initial_r::Float64, initial_eps_g::Float64)
    (initial_r <= 1.0) && error("initial_r must be greater than 1")
    (initial_eps_g <= 0.0) && error("eps_g must be greater than 0")
    return
  end

  function check_convergence(LBD::Float64, UBD::Float64, atol::Float64)
    if (abs(UBD - LBD) < atol)
      (verb == 2 || verb == 1) && println("Algorithm Converged")
      return true
    end
    return false
  end

  function llp2_update(llp2_obj_val::Float64, r::Float64, eps_g::Float64,
                       xbar::Vector{Float64}, pbar::Vector{Float64},
                       UBD::Float64, UBDg::Float64,
                       disc_set::Vector{Vector{Vector{Float64}}}, i::Int)
    xstar = fill(NaN,(length(xbar),))
    UBDtemp = UBDg
    eps_g_out = eps_g
    if (llp2_obj_val < 0.0)
      if (UBD <= UBDg)
        UBDtemp = UBD
        xstar[:] = xbar[:]
      end
      eps_g_out = eps_g/r
    else
      push!(disc_set[i], pbar)
    end
    return UBDtemp, xstar, eps_g_out
  end

  include("sip_explicit.jl")
  #include("sip_implicit.jl")
