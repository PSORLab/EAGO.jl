"""
    Variable_DR!(X::Array{Interval{Float64}},mult_lo::Vector{Float64},
                 mult_hi::Vector{Float64},LBD::Float64,UBD::Float64)

Uses duality multiplier from box constraints along with global upper bound to
tighten the size of the node, `Xt`, in place. Inputs:
* `X::Array{Interval{Float64}}`: Node over which to perform domain reduction
* `mult_lo::Vector{Float64}`: Multipliers for variable (low bound)
* `mult_hi::Vector{Float64}`: Multipliers for variable (high bound)
* `LBD::Float64` - Lower bound of problem
* `UBD::Float64` - Global upper bound for B&B algorithm
"""
function Variable_DR!(X::Vector{Interval{Float64}},
                      mult_lo::Vector{Float64},
                      mult_hi::Vector{Float64},
                      LBD::Float64,
                      UBD::Float64)

  lower_cut::Float64 = 0.0
  upper_cut::Float64 = 0.0
  if (LBD<=UBD)
    for i = 1:length(X)
      if (mult_lo[i] > 0.0)
        lower_cut = X[i].hi - (UBD-LBD)/mult_lo[i]
        if (lower_cut > X[i].lo)
          X[i] = Interval(lower_cut,X[i].hi)
        end
      elseif (mult_hi[i] > 0.0)
        upper_cut = X[i].lo + (UBD-LBD)/mult_lo[i]
        if (upper_cut < X[i].hi)
          X[i] = Interval(X[i].lo,upper_cut)
        end
      end
    end
  end
end
function Variable_DR!(X::Vector{MCInterval{Float64}},
                      mult_lo::Vector{Float64},
                      mult_hi::Vector{Float64},
                      LBD::Float64,
                      UBD::Float64)

  lower_cut::Float64 = 0.0
  upper_cut::Float64 = 0.0
  if (LBD<=UBD)
    for i = 1:length(X)
      if (mult_lo[i] > 0.0)
        lower_cut = X[i].hi - (UBD-LBD)/mult_lo[i]
        if (lower_cut > X[i].lo)
          X[i] = MCInterval(lower_cut,X[i].hi)
        end
      elseif (mult_hi[i] > 0.0)
        upper_cut = X[i].lo + (UBD-LBD)/mult_lo[i]
        if (upper_cut < X[i].hi)
          X[i] = MCInterval(X[i].lo,upper_cut)
        end
      end
    end
  end
end
"""
    Variable_DR_Imp!(X::Array{Interval{Float64}},mult_lo::Vector{Float64},
                     mult_hi::Vector{Float64},LBD::Float64,UBD::Float64,nx::Int64)

Uses duality multiplier from box constraints along with global upper bound to
tighten the size of the node, `Xt`, in place. Only tightens the X[(nx+1):(nx+np)]
intervals. Inputs:
* `X::Array{Interval{Float64}}`: Node over which to perform domain reduction
* `mult_lo::Vector{Float64}`: Multipliers for variable (low bound)
* `mult_hi::Vector{Float64}`: Multipliers for variable (high bound)
* `LBD::Float64` - Lower bound of problem
* `UBD::Float64` - Global upper bound for B&B algorithm
* `nx::Int64` - The number of state variables.
"""
function Variable_DR_Imp!(X::Vector{Interval{Float64}},
                      mult_lo::Vector{Float64},
                      mult_hi::Vector{Float64},
                      LBD::Float64,
                      UBD::Float64,
                      nx::Int64)

  lower_cut::Float64 = 0.0
  upper_cut::Float64 = 0.0
  if (LBD<=UBD)
    for i = (nx+1):length(X)
      if (mult_lo[i-nx] > 0.0)
        lower_cut = X[i].hi - (UBD-LBD)/mult_lo[i-nx]
        if (lower_cut > X[i].lo)
          X[i] = Interval(lower_cut,X[i].hi)
        end
      elseif (mult_hi[i-nx] > 0.0)
        upper_cut = X[i].lo + (UBD-LBD)/mult_lo[i-nx]
        if (upper_cut < X[i].hi)
          X[i] = Interval(X[i].lo,upper_cut)
        end
      end
    end
  end
end
function Variable_DR_Imp!(X::Vector{MCInterval{Float64}},
                      mult_lo::Vector{Float64},
                      mult_hi::Vector{Float64},
                      LBD::Float64,
                      UBD::Float64,
                      nx::Int64)

  lower_cut::Float64 = 0.0
  upper_cut::Float64 = 0.0
  if (LBD<=UBD)
    for i = (nx+1):length(X)
      if (mult_lo[i-nx] > 0.0)
        lower_cut = X[i].hi - (UBD-LBD)/mult_lo[i-nx]
        if (lower_cut > X[i].lo)
          X[i] = MCInterval(lower_cut,X[i].hi)
        end
      elseif (mult_hi[i-nx] > 0.0)
        upper_cut = X[i].lo + (UBD-LBD)/mult_lo[i-nx]
        if (upper_cut < X[i].hi)
          X[i] = MCInterval(X[i].lo,upper_cut)
        end
      end
    end
  end
end
