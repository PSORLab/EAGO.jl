# Solving Semi-Infinite Programs

## Using EAGO to Solve a Semi-Infinite Program

This example is also provided [here as a Jupyter Notebook](https://github.com/PSORLab/EAGO-notebooks/blob/master/notebooks/sip_explicit_solve.ipynb).

Semi-infinite programming remains an active area of research. In general, the solutions of semi-infinite programs (SIPs) with nonconvex semi-infinite constraints of the following form are extremely challenging:

```math
\begin{aligned}
f^{*} = & \min_{\mathbf x \in X} f(\mathbf x)\\
{\rm s.t.} \;\; & \max_{\mathbf p \in P} g(\mathbf x, \mathbf p) \leq 0\\
& \mathbf x \in X = \{ x \in \mathbb R : \mathbf x^{\mathbf L} \leq \mathbf x \leq \mathbf x^{\mathbf U} \}\\
& \mathbf p \in P = \{ p \in \mathbb R : \mathbf p^{\mathbf L} \leq \mathbf p \leq \mathbf p^{\mathbf U} \}
\end{aligned}
```

EAGO implements three different algorithms detailed in [[1](#References),[2](#References)] to determine a globally optimal solution to problems of the above form. This is accomplished using the [`sip_solve`](@ref) function which returns the optimal value, the solution, and a boolean feasibility flag. To illustrate the use of this function, a simple example is presented here which solves the problem:

```math
\begin{aligned}
f(\mathbf x) & = \frac{1}{3} x_{1}^{2} + x_{2}^{2} + \frac{x_{1}}{2}\\
g(\mathbf x, p) & = (1 - x_{1}^{2} p^{2})^{2} - x_{1} p^{2} - x_{2}^{2} + x_{2} \leq 0\\
& \mathbf x \in X = [-1000, 1000]^{2}\\
& p \in P = [0, 1]\\
\end{aligned}
```

```julia
using EAGO, JuMP

# Define semi-infinite program
f(x) = (1/3)*x[1]^2 + x[2]^2 + x[1]/2
gSIP(x, p) = (1.0 - (x[1]^2)*(p[1]^2))^2 - x[1]*p[1]^2 - x[2]^2 + x[2]

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
