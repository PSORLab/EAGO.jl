# Solving Semi-Infinite Programming

[Matthew Wilhelm](https://psor.uconn.edu/person/matthew-wilhelm/)  
Department of Chemical and Biomolecular Engineering, University of Connecticut

### Using EAGO to solve a SIP

Semi-infinite programming remains an active area of research. In general, the solution of semi-infinite programs with nonconvex semi-infinite constraints of the below formare extremely challenging:

$
\begin{align}
    f^=&\min_{\mathbf x}f(\mathbf{x}) \\
    &{\rm s.t.}\; \max_{\mathbf{p}}{g(\mathbf{x,p})}\le 0 \\
    & \mathbf{x} \in X = \{x \in \mathbb{R} : \mathbf{x^{L}} \leq \mathbf{x} \leq \mathbf{x^{U}}\} \\
    & \mathbf{p} \in P = \{p \in \mathbb{R} : \mathbf{p^{L}} \leq \mathbf{p} \leq \mathbf{p^{U}}\} \\
\end{align}
$

EAGO implements the SIPres of [1] to determine a globally optimal solution to problems of the above form. This accomplished using the `explicit_sip_solve` function which returns the optimal value, the solution, and a boolean feasibility value. To illustrate the functions use, a simple example is presented here which solves the below problem:

$
\begin{align}
    f(\mathbf{x}) &= (1/3)x_1^2 + x_2^2 + x_1/2 \\
    g(\mathbf{x},p) &= (1-(x_1^2)(p^2))^2 - x_1p^2 - x_2^2 + x_2 \leq 0 \\
     &{\; \qquad}\mathbf{x} \in X = [-1000, 1000]^2 \\
     &{\; \qquad}p \in P = [0.0, 1.0]  \label{ex:equalSIP}
\end{align}
$

```julia
using EAGO, JuMP

# Define semi-infinite program
f(x) = (1/3)*x[1]^2 + x[2]^2 + x[1]/2
gSIP(x,p) = (1.0 - (x[1]^2)*(p[1]^2))^2 - x[1]*p[1]^2 - x[2]^2 + x[2]

x_l = [-1000.0, -1000.0]
x_u = [1000.0, 1000.0]
p_l = [0.0]
p_u = [1.0]

# Create optimizer for use in solving SIP
opt = with_optimizer(EAGO.Optimizer, cut_max_iterations = 0,
                                         subgrad_tighten = false,
                                         absolute_tolerance = 1E-4,
                                         relative_tolerance = 1E-4,
                                         verbosity = 0)
m = JuMP.Model(opt)
sip_result = explicit_sip_solve(f, gSIP, x_l, x_u, p_l, p_u, m)

println("The global minimum of the semi-infinite program is between: $(sip_result.lower_bound) and $(sip_result.upper_bound).")
println("The global minimum is attained at: x = $(sip_result.xbar).")
println("Is the problem feasible? (sip_result.feasibility).")
```

### Using EAGO to solve a SIP w/ equality constraints

In general, the solution of semi-infinite programs with nonconvex semi-infinite equality in the form given below are extremely challenging. The general form of such problems is shown below:

$
\begin{align}\label{form:SIPImplicit}
    f^=&\min_{\mathbf y}f(\mathbf{y}) \\
    &{\rm s.t.}\; \max_{\mathbf{p,z}}{g(\mathbf{z,y,p})}\le 0 \\
    &{\; \qquad}{\rm s.t.}\; \mathbf{h(z,y,p) =  0} \\
    & \mathbf{y} \in Y = \{y \in \mathbb{R} : \mathbf{y^{L}} \leq \mathbf{y} \leq \mathbf{y^{U}}\} \\
    & \mathbf{p} \in P = \{p \in \mathbb{R} : \mathbf{p^{L}} \leq \mathbf{p} \leq \mathbf{p^{U}}\} \\
    & \mathbf{z} \in D_{x} \subset \mathbb{R}^{n_{x}}
\end{align}
$

EAGO implements the impSIPres algorithm of [1] to determine a globally optimal solution to problems of the above form subject to the stiputation that $\mathbf{h(z,y,p) =  0}$ uniquely determines $y$ as an implicit function of $z$ and $p$ on $Y \times P \times D_x$. This accomplished using the `implicit_sip_solve` function which returns the optimal value, the solution, and a boolean feasibility value. To illustrate the functions use, a simple example is presented here:

$
\begin{align}
    f(y) &= (y-3.5)^4 - 5(y-3.5)^3 -(y-3.5)^2 + 30(y-3.5) \\
    h(z,y,p) &= z - (y - y^3/6 + y^5/120)/\sqrt{z} - p = 0 \\
    g(z,y,p) &= z + cos(y-p/90) - p \leq 0 \\
     &{\; \qquad}y \in Y = [0.5,80]\\
     &{\; \qquad}p \in P = [80,120]  \label{ex:equalSIP}
\end{align}
$

```julia
using EAGO

# Variable bounds
xl = [0.5]; xu = [8.0]; yl = [68.8]; yu = [149.9]
pl = [80]; pu = [120];

# Objective and constraint functions are defined
f(y) =  (y[1]-3.5)^4 - 5(y[1]-3.5)^3 - (y[1]-3.5)^2 + 30(y[1]-3.5)
gSIP(z,y,p) = z[1] + cos(y[1]-p[1]/90) - p[1]
function h!(out,z,y,p)
    out[1] = z[1] - (y[1] - y[1]^3/6 + y[1]^5/120)/sqrt(z[1]) - p[1]
end

# The impSIPres routine is run
low_val, upp_val, sol, feas = implicit_sip_solve(f, h!, gSIP, xl, xu, yl, yu, pl, pu)

println("The global minimum of the semi-infinite program is between: $low_val and $upp_val.")
println("The global minimum is attained at: y = $val.")
```
