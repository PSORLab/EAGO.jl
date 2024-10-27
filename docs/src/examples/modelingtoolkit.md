# ModelingToolkit Example

This example is also provided [here as a Jupyter Notebook](https://github.com/PSORLab/EAGO-notebooks/blob/master/notebooks/mtk_system.ipynb).

### Using ModelingToolkit NonlinearSystem models with EAGO

This is a tutorial for exporting ModelingToolkit [[1](#References)] [`NonlinearSystem`](https://docs.sciml.ai/ModelingToolkit/stable/systems/NonlinearSystem/#ModelingToolkit.NonlinearSystem) models as standard Julia functions, which can then be used as EAGO-compatible equality constraints in JuMP models.

## Defining Nonlinear System

The system of interest is derived from an example originally presented by [[2](#References)] that involves a continuous stirred-tank reactor (CSTR) and separator train (with recycle) for the chlorination of benzene with the following reactions taking place:

```math
\begin{array}{rcl}
\text{C}_{6} \text{H}_{6} + \text{Cl}_{2} &\rightarrow& \text{C}_{6} \text{H}_{5} \text{Cl} + \text{HCl}\\
\text{C}_{6} \text{H}_{5} \text{Cl} + \text{Cl}_{2} &\rightarrow& \text{C}_{6} \text{H}_{4} \text{Cl}_{2} + \text{HCl}
\end{array}
```

where the rate constants ``k_{1}`` and ``k_{2}`` ``[h^{-1}]`` are known and the reactor volume ``V`` ``[m^{3}]`` and feed flow rate ``F_{1}`` ``[kmol/h]`` are considered free design variables. The CSTR is followed by a separation train for product purification and reactant recycle.

![mtk_pfd](assets/mtk_pfd.png)

For simplicity, the reactions will be first-order (no dependence on ``\text{Cl}_{2}``) with respect to benzene (A) and chlorobenzene (B). The molar volumes of each species are: ``V_{A} = 8.937 \times 10^{-2} \; m^{3}/kmol``, ``V_{B} = 1.018 \times 10^{-1} \;  m^{3}/kmol``, and ``V_{C} = 1.13 \times 10^{-1} \; m^{3}/kmol`` with dichlorobenzene as C. The feed is considered to be pure A (ignoring ``\text{Cl}_{2}``). The rate constants are ``k_{1} = 0.40 \; h^{-1}`` and ``k_{2} = 0.055 \; h^{-1}``. The steady-state model equations are:

```math
\begin{array}{rclr}
    F_{1} + F_{7} &=& F_{2} & (\text{mixer})\\	
	y_{3,A} F_{3} &=& F_{2} - r_{1} V & (\text{reactor})\\
	y_{3,B} F_{3} &=&  (r_{1} - r_{2}) V & (\text{reactor})\\
	y_{3,C} F_{3} &=&  r_{2} V & (\text{reactor})\\
	F_{3} &=& F_{4} + F_{7} & (\text{separator 1})\\
	y_{3,B} F_{3} &=& y_{4,B} F_{4} & (\text{separator 1})\\
	y_{3,C} F_{3} &=& y_{4,C} F_{4} & (\text{separator 1})\\
	F_{4} &=& F_{5} + F_{6} & (\text{separator 2})\\
	y_{4,B} F_{4} &=& F_{5} & (\text{separator 2})\\
	y_{3,A} + y_{3,B} + y_{3,C} &=& 1 & (\text{relationship})\\
	y_{4,B} + y_{4,C} &=& 1 & (\text{relationship})
\end{array}
```

where ``y_{i,j}`` is the mole fraction of ``j \in \{A,B,C\}`` in Stream ``i \in \{3,4\}`` (with ``y_{4,A} = 0``) and ``r_{1}`` and ``r_{2}`` are the reaction rates ``[kmol/(m^{3} h)]`` defined as:

```math
\begin{array}{rcl}
r_{1} &=& k_{1} y_{3,A}/(y_{3,A} V_{A} + y_{3,B} V_{B} + y_{3,C} V_{C})\\
r_{2} &=& k_{2} y_{3,B}/(y_{3,A} V_{A} + y_{3,B} V_{B} + y_{3,C} V_{C})
\end{array}
```

with the variable vector ``\mathbf{x} \in \mathbb{R}^{12}``:

```math
\mathbf{x} = (V, F_{1}, F_{2}, y_{3,A}, y_{3,B}, y_{3,C}, F_{3}, y_{4,B}, y_{4,C}, F_{4}, F_{6}, F_{7}).
```

### Performance Specifications

We must produce at least 25 ``kmol/h`` of B (``F_{5} \ge 25``) and cannot have a reactor larger than 10 ``m^{3}`` due to space limitations. From a laboratory study, it was found that the residence time ``\tau`` ``[h]`` for the reactor must be at least 475 seconds. Residence time can be defined as:

```math
\tau = \frac{V}{F_{3} (y_{3,A} V_{A} + y_{3,B} V_{B} + y_{3,C} V_{C})}.
```

### Design Problem

We will consider an economic objective for our optimization problem as the total annualized costs of the unit operations. The total annualized cost of the CSTR is given by

```math
f_{CSTR} = (25764 + 8178V)/2.5
```

and the total annualized cost of the separator train is given by:

```math
\begin{array}{rcl}
s_{1}^{cap} &=& 132718 + F_{3} (369 y_{3,A} - 1113.9 y_{3,B})\\
s_{2}^{cap} &=& 25000 + F_{4} (6984.5 y_{4,B} - 3869.53 y_{4,C}^{2})\\
s_{1}^{op} &=& F_{3} (3 + 36.11 y_{3,A} + 7.71 y_{3,B}) (26.32 \times 10^{-3})\\
s_{2}^{op} &=& F_{4} (26.21 + 29.45 y_{4,B}) (26.32 \times 10^{-3})\\
f_{sep} &=& (s_{1}^{cap} + s_{2}^{cap})/2.5 + 0.52 (s_{1}^{op} + s_{2}^{op}).
\end{array}
```

The total annualized cost of the complete system is then given by

```math
f_{total} = f_{CSTR} + f_{sep}.
```

## Modeling and Simplifying Using ModelingToolkit

```julia
using EAGO, GLPK, JuMP, ModelingToolkit
```

First, we define the state variables with `@variables` and model parameters with `@parameters`. Note that the design variables ``V`` and ``F_{1}`` are considered parameters so we have a square system (10 equations, 10 unknowns).

```julia
# State variables
ModelingToolkit.@variables F₂ y_3A y_3B y_3C F₃ y_4B y_4C F₄ F₆ F₇

# Model parameters and design variables
ModelingToolkit.@parameters k₁ k₂ V_A V_B V_C V F₁
```

Next, we can write symbolic expressions for ``r_{1}``, ``r_{2}``, and ``F_{5}`` to help simplify our model equations.

```julia
# Symbolic expressions for reaction rates and F₅
r₁ = (k₁*y_3A)/(y_3A*V_A + y_3B*V_B + y_3C*V_C)
r₂ = (k₂*y_3B)/(y_3A*V_A + y_3B*V_B + y_3C*V_C)
F₅ = (y_4B*F₄)
```

The variables, parameters, and steady-state model equations are then used to construct the `NonlinearSystem`.

```julia
# Steady-state model equations
eqs = [
    F₁ + F₇ ~ F₂
    (y_3A*F₃) ~ F₂ - r₁*V
    (y_3B*F₃) ~ (r₁ - r₂)*V
    (y_3C*F₃) ~ r₂*V
    F₃ ~ F₄ + F₇
    (y_3B*F₃) ~ (y_4B*F₄)
    (y_3C*F₃) ~ (y_4C*F₄)
    F₄ ~ F₅ + F₆
    y_3A + y_3B + y_3C ~ 1
    y_4B + y_4C ~ 1
]

# Variables and parameters
vars = [F₂, y_3A, y_3B, y_3C, F₃, y_4B, y_4C, F₄, F₆, F₇]
pars = [k₁, k₂, V_A, V_B, V_C, V, F₁]

# Building and simplifying model
@mtkbuild ns = NonlinearSystem(eqs, vars, pars)
```

`@mtkbuild` simplifies the model down to 4 equations and 4 unknowns, significantly reducing the dimensionality of this system. The variables that were eliminated during the simplification process are called "observed variables," and their expressions can be obtained using `observed(::NonlinearSystem)`.

## Converting to a Standard Function

Now that we've simplified our model, we need to covert the symbolic model equations into standard Julia functions in order to use them with EAGO. This can be accomplished through the use of `Symbolics.build_function`, which generates a numerically-usable function from a symbolic expression. Since `Symbolics.jl` is already built into `ModelingToolkit.jl`, no additional packages need to be installed.

To simplify the use of `build_function` for our purposes, we first create the function `ns_to_function()` that takes in a `NonlinearSystem`, a vector of parameters, and a vector of function inputs, and returns a vector of standard Julia functions corresponding to the simplified model equations in the nonlinear system.

```julia
# Converts symbolic model equations into standard Julia functions
function ns_to_function(ns::NonlinearSystem, params::Vector{Pair{Num,Float64}}, vars::Vector{Num})
    function_holder = []
    # For each model equation, substitute parameter values and create function
    for eqn in full_equations(ns)
        expr = substitute(eqn.rhs, Dict{Num,Any}(params))
        new_func = build_function(expr, vars..., expression = Val{false})
        push!(function_holder, new_func)
    end
    return function_holder
end
```

To use this function, we first assign values to our known parameters, define our function inputs, and then call `ns_to_function()` to obtain our desired functions.

```julia
# Assign parameter values
params = [
    k₁ => 0.40,
    k₂ => 0.055,
    V_A => 8.937e-2,
    V_B => 1.018e-1,
    V_C => 1.130e-1
]

# Function inputs (design and state variables)
x = [V, F₁, y_3B, y_3C, F₃, y_4C]

# Create functions for each model equation as a function of design and state variables
f = ns_to_function(ns, params, x)
```

We now have standard Julia functions for each model equation which can be used as equality constraints in our JuMP [[3](#References)] model. 

We'll also create functions for our inequality constraints (``F_{5}``, ``\tau``) and objective function (``f_{total}``). Since these functions are not expressed in terms of the two design variables (``V`` and ``F_{1}``) and the four state variables (``y_{3,B}``, ``y_{3,C}``, ``F_{3}``, ``y_{4,C}``), we'll create a function `expr_to_function()` to substitute the `NonlinearSystem` observed variables into our input expression before calling `build_function`.

All we have to do now is write the symbolic expressions for ``F_{5}``, ``\tau``, and ``f_{total}``, then call `expr_to_function()` for each expression.

```julia
# Converts symbolic expression into standard Julia function, substituting parameter and observed variable expressions
function expr_to_function(expr::Num, ns::NonlinearSystem, params::Vector{Pair{Num,Float64}}, vars::Vector{Num})
    # Creates a dictionary of parameter and observed variable substitutions
    substitute_dict = Dict{Num,Any}(params)
    for eqn in observed(ns)
        substitute_dict[eqn.lhs] = eqn.rhs
    end
    # Makes substitutions until no substitutions can be made
    while ~isempty(intersect(string.(get_variables(expr)), string.(keys(substitute_dict))))
        expr = substitute(expr, substitute_dict)
    end
    return build_function(expr, vars..., expression = Val{false})
end

# Symbolic expression for F₅
exprF5 = y_4B*F₄

# Symbolic expression for τ
exprTau = V/(F₃*(y_3A*V_A + y_3B*V_B + y_3C*V_C))

# Symbolic expressions for CSTR and separator costs
f_CSTR = (25764 + 8178*V)/2.5
s1cap = 132718 + F₃*(369*y_3A - 1113.9*y_3B)
s2cap = 25000 + F₄*(6984.5*y_4B - 3869.53*y_4C^2)
s1op = F₃*(3+36.11*y_3A + 7.71*y_3B)*26.32e-3
s2op = F₄*(26.21 + 29.45*y_4B)*26.32e-3;
f_Sep = (s1cap+s2cap)/2.5 + 0.52*(s1op+s2op)

# Symbolic expression for total cost (objective function)
exprTotCost = f_CSTR + f_Sep

# Create functions for each expression as a function of design and state variables
F5 = expr_to_function(exprF5, ns, params, x)
Tau = expr_to_function(exprTau, ns, params, x)
TotCost = expr_to_function(exprTotCost, ns, params, x)
```

## Solving the Optimization Problem with EAGO

To formulate and solve the design problem, we'll begin by creating a JuMP model. We'll also use GLPK as a subsolver because it greatly accelerates convergence for this particular system.

```julia
# Using GLPK subsolver for EAGO
factory = () -> EAGO.Optimizer(SubSolvers(; r = GLPK.Optimizer()))
m = Model(optimizer_with_attributes(factory,
                                    "branch_cvx_factor" => 0.75,
                                    "output_iterations" => 200))
```

Next, we need to define bounds for our design variables (``V``, ``F_{1}``) and state variables (``y_{3,B}``, ``y_{3,C}``, ``F_{3}``, ``y_{4,C}``). This gives us a total of 6 decision variables in our optimization problem. 

```julia
# Define variable bounds
xL = [5.0, 25.0, 0.1, 0.001, 75.0, 0.01]
xU = [10.0, 50.0, 0.5, 0.1, 125.0, 0.1]
# x = (V, F₁, y_3B, y_3C, F₃, y_4C)
@variable(m, xL[i] <= x[i=1:6] <= xU[i])
```

!!! note

    It's important to note that defining "good" bounds is especially important for deterministic global solvers like EAGO because of the set arithmetic and interval and convex analyses required to guarantee global optimality. Based on simple linear constraints and understanding what our variables physically represent, we can easily eliminate erroneous regions from the search space and drastically reduce computational cost. For example, we know the reactor volume must be less than or equal to 10 ``m^{3}``, but we also know it should be at least a few cubic meters large, so we bound ``V`` between 5 and 10. Similarly, we know we must produce at least 25 ``kmol/h`` of B, but this is only possible if there are at least 25 ``kmol/h`` of A entering the system since we are operating under steady-state conditions with one influent and two effluent streams.

Once proper bounds are defined, we continue by registering our user-defined functions and defining our nonlinear constraints (model equations and performance specifications) and objective (minimizing total annualized costs).

```julia
# Register user-defined functions
JuMP.register(m, :f1, 6, f[1], autodiff=true)
JuMP.register(m, :f2, 6, f[2], autodiff=true)
JuMP.register(m, :f3, 6, f[3], autodiff=true)
JuMP.register(m, :f4, 6, f[4], autodiff=true)
JuMP.register(m, :F5, 6, F5, autodiff=true)
JuMP.register(m, :Tau, 6, Tau, autodiff=true)
JuMP.register(m, :TotCost, 6, TotCost, autodiff=true)

# Define model equation constraints
@NLconstraint(m, e1, f1(x...) == 0)
@NLconstraint(m, e2, f2(x...) == 0)
@NLconstraint(m, e3, f3(x...) == 0)
@NLconstraint(m, e4, f4(x...) == 0)

# Define performance specification constraints
@NLconstraint(m, c1, F5(x...) ≥ 25.0)
@NLconstraint(m, c2, Tau(x...) ≥ 475.0/3600.0)

# Define objective to minimize total annualized costs
@NLobjective(m, Min, TotCost(x...))
```

We're now ready to solve the optimal design problem. We'll `@time` the solver so we can compare convergence times for reduced-space vs full-space formulations later.

```julia
# Solve optimization problem
@time JuMP.optimize!(m)
```

We've identified a global optimal solution ``\mathbf{x}^* = (8.4594, 26.3167, 0.2631, 0.01386, 95.0215, 0.05003)``. Keep in mind that we can readily calculate values for the observed variables since we have their explicit expressions written in terms of these decision variables.

## Full-Space Formulation

To compare and showcase the full-space optimization formulation for this system, we'll quickly formulate and solve the design problem with all 10 model equations using the full variable vector

```math
\mathbf{x} = (V, F_{1}, F_{2}, y_{3,A}, y_{3,B}, y_{3,C}, F_{3}, y_{4,B}, y_{4,C},F_{4}, F_{6}, F_{7}).
```

We can generate functions for each of the original model equations by first initializing the model using `@named` instead of `@mtkbuild` to exclude the structural simplification step.

```julia
# Initialize model without structural simplification
@named fns = NonlinearSystem(eqs, vars, pars)
```

Then, we redefine our function inputs and generate our functions in the same way as before.

```julia
# Full variable vector
x = [V, F₁, F₂, y_3A, y_3B, y_3C, F₃, y_4B, y_4C, F₄, F₆, F₇]

# Create functions for each model equation as a function of the full variable vector
ff = ns_to_function(fns, params, x)

# Create functions for each expression as a function of the full variable vector
fF5 = expr_to_function(exprF5, fns, params, x)
fTau = expr_to_function(exprTau, fns, params, x)
fTotCost = expr_to_function(exprTotCost, fns, params, x)
```

We'll define a new JuMP model with the same settings, define bounds for each variable, register functions, define constraints and objective, then solve the optimization problem.

```julia
n = Model(optimizer_with_attributes(factory,
                                    "branch_cvx_factor" => 0.75,
                                    "output_iterations" => 200))

# Define variable bounds
# New variables:   v    v                      v           v    v     v
xL = [ 5.0, 25.0, 75.0, 0.5, 0.1, 0.001, 75.0, 0.9, 0.01, 25.0, 0.0, 50.0]
xU = [10.0, 50.0, 125.0, 1.0, 0.5, 0.1, 125.0, 1.0, 0.1, 50.0, 10.0, 100.0]
# x = (V, F₁, F₂, y_3A, y_3B, y_3C, F₃, y_4B, y_4C, F₄, F₆, F₇)
@variable(n, xL[i] <= x[i=1:12] <= xU[i])

# Register user-defined functions
JuMP.register(n, :f1, 12, ff[1], autodiff=true)
JuMP.register(n, :f2, 12, ff[2], autodiff=true)
JuMP.register(n, :f3, 12, ff[3], autodiff=true)
JuMP.register(n, :f4, 12, ff[4], autodiff=true)
JuMP.register(n, :f5, 12, ff[5], autodiff=true)
JuMP.register(n, :f6, 12, ff[6], autodiff=true)
JuMP.register(n, :f7, 12, ff[7], autodiff=true)
JuMP.register(n, :f8, 12, ff[8], autodiff=true)
JuMP.register(n, :f9, 12, ff[9], autodiff=true)
JuMP.register(n, :f10, 12, ff[10], autodiff=true)
JuMP.register(n, :F5, 12, fF5, autodiff=true)
JuMP.register(n, :Tau, 12, fTau, autodiff=true)
JuMP.register(n, :TotCost, 12, fTotCost, autodiff=true)

# Define model equation constraints
@NLconstraint(n, e1, f1(x...) == 0)
@NLconstraint(n, e2, f2(x...) == 0)
@NLconstraint(n, e3, f3(x...) == 0)
@NLconstraint(n, e4, f4(x...) == 0)
@NLconstraint(n, e5, f5(x...) == 0)
@NLconstraint(n, e6, f6(x...) == 0)
@NLconstraint(n, e7, f7(x...) == 0)
@NLconstraint(n, e8, f8(x...) == 0)
@NLconstraint(n, e9, f9(x...) == 0)
@NLconstraint(n, e10, f10(x...) == 0)

# Define performance specification constraints
@NLconstraint(n, c1, F5(x...) ≥ 25.0)
@NLconstraint(n, c2, Tau(x...) ≥ 475.0/3600.0)

# Define objective to minimize total annualized costs
@NLobjective(n, Min, TotCost(x...))

# Solve optimization problem
@time JuMP.optimize!(n)
```

It takes three times as many iterations and significantly more memory allocations for EAGO to converge with the full-space formulation, demonstrating the advantages of reduced-space formulations. Because global solvers suffer from the curse of dimensionality, reducing the number of decision variables through model simplifications can improve performance and therefore broaden the scope of problems for which global solvers can be applied.

This tutorial introduces user-defined functions that provide an easy way to use ModelingToolkit's intuitive modeling interface and model simplification features to create reduced-space formulations in a way that is compatible with EAGO.

## References

1. Yingbo Ma and Shashi Gowda and Ranjan Anantharaman and Chris Laughman and Viral Shah and Chris Rackauckas, ModelingToolkit: A Composable Graph Transformation System For Equation-Based Modeling, *arXiv*, (2021).
2. A. C. Kokossis, C. A. Floudas. Synthesis of isothermal reactor—separator—recycle systems, *Chemical Engineering Science*, 46:5-6 (1991), pp. 1361-1383.
3. Iain Dunning and Joey Huchette and Miles Lubin. JuMP: A Modeling Language for Mathematical Optimization, *SIAM Review*, 59 (2017), pp. 295-320.