function variable_duality_based_tightening!(x::NodeBB,
                        mult_lo::Vector{Float64},
                        mult_hi::Vector{Float64},
                        LBD::Float64,
                        UBD::Float64)
  nx::Int = length(x)
  lower_cut::Float64 = 0.0
  upper_cut::Float64 = 0.0
  if (LBD <= UBD)
    for i = 1:nx
      if (mult_lo[i] > 0.0)
        lower_cut = x.upper_variable_bounds[i] - (UBD - LBD)/mult_lo[i]
        (lower_cut > x.lower_variable_bounds[i]) && (x.lower_variable_bounds[i] = lower_cut)
      elseif (mult_hi[i] > 0.0)
        upper_cut = x.lower_variable_bounds[i] + (UBD - LBD)/mult_lo[i]
        (upper_cut < x.upper_variable_bounds[i]) && (x.upper_variable_bounds[i] = upper_cut)
      end
    end
  end
end
