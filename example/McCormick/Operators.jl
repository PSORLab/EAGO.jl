workspace()

# x*(x-3)*sin(x*5) on 1 to 4

using EAGO

# create seed gradient
xSMCg = SMCg{1,Interval{Float64},Float64}(2.0,                   # Concave relaxation
                                          2.0,                   # Convex relaxation
                                          seed_g(Float64,1,1),   # Create concave subgradient
                                          seed_g(Float64,1,1),   # Create convex subgradient
                                          Interval(1.0,4.0),     # Set interval bounds
                                          false)

relaxed_f = xSMCg*(xSMCg-3)*sin(xSMCg*5)
