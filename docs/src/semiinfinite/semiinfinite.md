# Solving Semi-Infinite Programs

## Using EAGO to Solve a Semi-Infinite Program (SIP)

Semi-infinite programming remains an active area of research. In general, the solutions of SIPs with nonconvex semi-infinite constraints of the following form are extremely challenging:

![SipProbForm](SIPProbFormulation.png)

EAGO implements three different algorithm detailed in [[1](#references),[2](#references)] to determine a globally optimal solution to problems of the above form. This is accomplished using the `sip_solve` function which returns the optimal value, the solution, and a boolean feasibility flag. To illustrate the use of this function, a simple example is presented here which solves the problem:

![SipForm](SIPformulation.png)

```julia
using EAGO, JuMP

# Define semi-infinite program
f(x) = (1/3)*x[1]^2 + x[2]^2 + x[1]/2
gSIP(x,p) = (1.0 - (x[1]^2)*(p[1]^2))^2 - x[1]*p[1]^2 - x[2]^2 + x[2]

x_l = Float64[-1000.0, -1000.0]
x_u = Float64[1000.0, 1000.0]
p_l = Float64[0.0]
p_u = Float64[1.0]

sip_result = sip_solve(SIPRes(), x_l, x_u, p_l, p_u, f, Any[gSIP], abs_tolerance = 1E-3)

println("The global minimum of the semi-infinite program is between: $(sip_result.lower_bound) and $(sip_result.upper_bound).")
println("The global minimum is attained at: x = $(sip_result.xsol).")
println("Is the problem feasible? $(sip_result.feasibility).")
```

## Semi-Infinite Solver

```@docs
    SIPProblem
    SIPResult
    SIPRes
    SIPResRev
    SIPHybrid
    get_sip_optimizer
    build_model
    sip_llp!
    sip_bnd!
    sip_res!
    sip_solve
```

## References

1. **Mitsos A (2009).** Global optimization of semi-infinite programs via restriction of the right-hand side. *Optimization*, 60(10-11):1291-1308.
2. **Djelassi, Hatim, and Alexander Mitsos.** A hybrid discretization algorithm with guaranteed feasibility for the global solution of semi-infinite programs. *Journal of Global Optimization*, 68.2 (2017): 227-253 should be used.
