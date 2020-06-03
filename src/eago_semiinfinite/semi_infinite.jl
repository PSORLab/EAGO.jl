# Copyright (c) 2018: Matthew Wilhelm & Matthew Stuber.
# This work is licensed under the Creative Commons Attribution-NonCommercial-
# ShareAlike 4.0 International License. To view a copy of this license, visit
# http://creativecommons.org/licenses/by-nc-sa/4.0/ or send a letter to Creative
# Commons, PO Box 1866, Mountain View, CA 94042, USA.
#############################################################################
# EAGO
# A development environment for robust and global optimization
# See https://github.com/PSORLab/EAGO.jl
#############################################################################
# TODO
#############################################################################

"""
      SIPResult

Structure storing the results of the SIPres algorithm.
"""
  mutable struct SIPResult
      iteration_number::Int64
      upper_bound::Float64
      lower_bound::Float64
      feasibility::Bool
      xsol::Vector{Float64}
      psol::Vector{Float64}
      solution_time::Float64
  end
  SIPResult() = SIPResult(1, Inf, -Inf, true, Float64[], Float64[], 0.0)

"""
      SIPProblem

Structure storing problem information for the solution routine.
"""
mutable struct SIPProblem
    x_l::Vector{Float64}
    x_u::Vector{Float64}
    p_l::Vector{Float64}
    p_u::Vector{Float64}

    np::Int64
    nSIP::Int64
    nx::Int64
    sense::Symbol

    init_lower_disc::Vector{Vector{Vector{Float64}}}
    init_upper_disc::Vector{Vector{Vector{Float64}}}

    absolute_tolerance::Float64
    constraint_tolerance::Float64
    iteration_limit::Int64
    initial_eps_g::Float64
    initial_r::Float64

    return_hist::Bool
    header_interval::Int64
    print_interval::Int64
    verbosity::Int64

    optimizer
    kwargs
  end
  function SIPProblem(x_l::Vector{Float64}, x_u::Vector{Float64},
                      p_l::Vector{Float64}, p_u::Vector{Float64},
                      gSIP, optimizer, kwargs)

      np = length(p_l)
      nx = length(x_l)

      absolute_tolerance = haskey(kwargs, :sip_absolute_tolerance) ? kwargs[:sip_absolute_tolerance] : 1E-3
      constraint_tolerance = haskey(kwargs, :sip_constraint_tolerance) ? kwargs[:sip_constraint_tolerance] : 1E-3
      iteration_limit = haskey(kwargs, :sip_iteration_limit) ? kwargs[:sip_iteration_limit] : 100
      initial_eps_g = haskey(kwargs, :sip_initial_eps_g) ? kwargs[:sip_initial_eps_g] : 1.0
      initial_r = haskey(kwargs, :sip_initial_r) ? kwargs[:sip_initial_r] : 2.0

      return_hist = haskey(kwargs, :sip_return_hist) ? kwargs[:sip_return_hist] : false
      header_interval = haskey(kwargs, :sip_header_interval) ? kwargs[:sip_header_interval] : 20
      print_interval = haskey(kwargs, :sip_print_interval) ? kwargs[:sip_print_interval] : 1
      verbosity = haskey(kwargs, :sip_verbosity) ? kwargs[:sip_verbosity] : 1

      sense = haskey(kwargs, :sip_sense) ? kwargs[:sip_sense] : :min
      init_lower_disc = haskey(kwargs, :sip_init_lower_disc) ? kwargs[:sip_init_lower_disc] : Vector{Vector{Float64}}[]
      init_upper_disc = haskey(kwargs, :sip_init_upper_disc) ? kwargs[:sip_init_upper_disc] : Vector{Vector{Float64}}[]

      opt_dict = Dict{Symbol,Any}()
      for key in keys(kwargs)
        string_key = String(key)
        if string_key[1:3] !== "sip"
          opt_dict[key] = kwargs[key]
        end
      end

      nSIP = length(gSIP)

      SIPProblem(x_l, x_u, p_l, p_u, np, nSIP, nx, sense, init_lower_disc,
                 init_upper_disc, absolute_tolerance, constraint_tolerance,
                 iteration_limit,
                 initial_eps_g, initial_r, return_hist, header_interval,
                 print_interval, verbosity, optimizer, opt_dict)
  end

  function print_int!(verbosity::Int64, hdr_intv::Int64, prt_intv::Int64,
                      k_int::Int64, lbd::Float64, ubd::Float64,
                      eps::Float64, r::Float64, ismin::Bool)

    if (verbosity == 1 || verbosity == 2)

      # prints header line every hdr_intv times
      if (mod(k_int, hdr_intv) == 0 || k_int == 1)
        println("| Iteration | Lower Bound | Upper Bound |   eps   |   r   |  Gap  |  Ratio  |")
      end

      # prints iteration summary every prnt_intv times
      if mod(k_int, prt_intv) == 0

        print_str = "| "

        max_len = 15
        temp_str = string(k_int)
        len_str = length(temp_str)
        print_str *= (" "^(max_len - len_str))*temp_str*" | "

        max_len = 15
        lower_adj = ismin ? lbd : ubd
        temp_str = string(round(lower_adj, digits = 6))
        len_str = length(temp_str)
        print_str *= (" "^(max_len - len_str))*temp_str*" | "

        max_len = 15
        upper_adj = ismin ? ubd : lbd
        temp_str = string(round(upper_adj, digits = 6))
        len_str = length(temp_str)
        print_str *= (" "^(max_len - len_str))*temp_str*" | "

        max_len = 15
        temp_str = string(round(eps, digits = 6))
        len_str = length(temp_str)
        print_str *= (" "^(max_len - len_str))*temp_str*" | "

        max_len =15
        temp_str = string(round(r, digits = 5))
        len_str = length(temp_str)
        print_str *= (" "^(max_len - len_str))*temp_str*" | "

        max_len = 15
        temp_str = string(round((ubd-lbd), digits = 5))
        len_str = length(temp_str)
        print_str *= (" "^(max_len - len_str))*temp_str*" | "

        max_len = 15
        temp_str = string(round((ubd-lbd)/abs(ubd), digits = 6))
        len_str = length(temp_str)
        print_str *= (" "^(max_len - len_str))*temp_str*" |"

        println(print_str)
      end
    end
    return
  end

  function print_summary!(verb::Int64, val::Float64, x::Vector{Float64},
                          feas::Bool, desc::String)
    (verb == 2 || verb == 1) && println("solved $desc: ", val, " ", x, " ", feas)
    return
  end

  function check_inputs!(initial_r::Float64, initial_eps_g::Float64)
    (initial_r <= 1.0) && error("initial_r must be greater than 1")
    (initial_eps_g <= 0.0) && error("eps_g must be greater than 0")
    return
  end

  function check_convergence(LBD::Float64, UBD::Float64, atol::Float64, verb::Int64)
    if (abs(UBD - LBD) < atol)
      (verb == 2 || verb == 1) && println("Algorithm Converged")
      return true
    end
    return false
  end

  struct SIPCallback
    f
    gSIP
  end

  include("sip_explicit.jl")
