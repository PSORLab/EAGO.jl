# Advanced-Use Example 1

This example is also provided [here as a Jupyter Notebook](https://github.com/PSORLab/EAGO-notebooks/blob/master/notebooks/custom_quasiconvex.ipynb).

### Customizing EAGO to Solve a Quasiconvex Problem

In this example, we'll adapt EAGO to implement the bisection-based algorithm used to solve a quasiconvex optimization problem presented in [[1](#References)]:

```math
\begin{aligned}
f^{*} = & \min_{\mathbf y \in Y} f(\mathbf y) \\
{\rm s.t.} \; \; & \sum_{i = 1}^{5} i \cdot y_{i} - 5 = 0 \\
& \sum_{i = 1}^{5} y_{i}^{2} - 0.5 \pi \leq 0 \\
& -\bigg(\frac{1}{2} y_{1}^{2} + \frac{1}{2} y_{2}^{2} + 2 y_{1} y_{2} + 4 y_{1} y_{3} + 2 y_{2} y_{3} \bigg) \leq 0 \\
& -y_{1}^{2} - 6 y_{1} y_{2} - 2 y_{2}^{2} + \cos (y_{1}) + \pi \leq 0 \\
& Y = [0, 5]^{5}
\end{aligned}
```

where

```math
\begin{aligned}
f(\mathbf y) = -\frac{\ln ((5 + y_{1})^{2} + \sum_{i = 1}^{5} y_{i})}{1 + \sum_{i = 1}^{5} y_{i}^{2}}.
\end{aligned}
```

Interval analysis shows that the objective value is bounded by the interval ``F`` such that ``f^{*} \in F = [f^{L}, f^{U}] = [-5, 0]``. Introducing an auxiliary variable ``t \in T = F`` allows the problem to be formulated as:

```math
\begin{aligned}
t^{*} = & \min_{\mathbf y \in Y, t \in T} t \\
{\rm s.t.} \; \; & (24) - (27) \\
& f(\mathbf y) - t \leq 0 \\
& Y = [0,5]^{2}, \; \; T = [-5,0].
\end{aligned}
```

Let ``\phi_{\tau}(\mathbf y) = f(\mathbf y) - \tau`` such that ``\tau = (t^{L} + t^{U})/2``. We solve for ``\mathbf y`` subject to constraints ``(24) - (27)`` where ``\phi_{\tau} (\mathbf y) \leq 0``. If this is feasible, ``t^{*} \in [t^{L},\tau]``, else ``t^{*} \in [τ, t^{U}]``. The interval containing ``t^{*}`` is kept and the other is fathomed. This manner of bisection is repeated until an interval containing a feasible solution with a width of at most ``\epsilon`` is located [[2](#References)].

## Customizing EAGO's Script

First, the preprocessing step, upper problem, and postprocessing routines are short-circuited as only a single optimization problem needs to be solved at each iteration.

```julia
using MathOptInterface, EAGO, JuMP
import EAGO: Optimizer, GlobalOptimizer

struct QuasiConvex <: EAGO.ExtensionType end
import EAGO: preprocess!, upper_problem!, postprocess!
function EAGO.preprocess!(t::QuasiConvex, x::GlobalOptimizer)
    x._preprocess_feasibility = true
end
function EAGO.upper_problem!(t::QuasiConvex, x::GlobalOptimizer)
    x._upper_feasibility = true
end
function EAGO.postprocess!(t::QuasiConvex, x::GlobalOptimizer)
    x._postprocess_feasibility = true
end
```

Next, we specify that only an absolute tolerance should be checked for convergence and termination.

```julia
import EAGO: convergence_check, termination_check
function EAGO.convergence_check(t::QuasiConvex, x::GlobalOptimizer)
    gap = (x._upper_objective_value - x._lower_objective_value)
    return (gap <= x._parameters.absolute_tolerance)
end
function EAGO.termination_check(t::QuasiConvex, x::GlobalOptimizer)
    flag = EAGO.convergence_check(t, x)
    if flag
        x._end_state = EAGO.GS_OPTIMAL
        x._termination_status_code = MathOptInterface.OPTIMAL
        x._result_status_code = MathOptInterface.FEASIBLE_POINT
    end
    return flag
end
```

We then indicate that only the sixth variable, representing ``t``, should be branched on. Since we will apply our knowledge about which ``t^{*}`` should be kept in the lower problem definition, we also short-circuit the [`EAGO.repeat_check`](@ref) function here to tell EAGO not to branch this node, but instead to repeatedly evaluate it.

```julia
import EAGO: repeat_check
branch_variable = [i == 6 for i=1:6]
EAGO.repeat_check(t::QuasiConvex, x::GlobalOptimizer) = true
```

In the lower problem, we then specify that the problem is to be solved locally for a fixed ``t`` value. The objective value is then updated and the problem is contracted in order to discard the region which is known to not contain the optimal value.

```julia
import EAGO: lower_problem!
function EAGO.lower_problem!(t::QuasiConvex, x::GlobalOptimizer)
    y = x._current_node
    indx = x._sol_to_branch_map[6]
    lower = y.lower_variable_bounds[indx]
    upper = y.upper_variable_bounds[indx]
    midy = (lower + upper)/2.0
    y.lower_variable_bounds[indx] = midy
    y.upper_variable_bounds[indx] = midy
    EAGO.solve_local_nlp!(x)
    feas = x._upper_feasibility
    y.lower_variable_bounds[indx] = feas ? lower : midy
    y.upper_variable_bounds[indx] = feas ? midy : upper
    x._lower_objective_value = y.lower_variable_bounds[indx]
    x._upper_objective_value = y.upper_variable_bounds[indx]
    x._lower_feasibility = true
    return
end
```

We now define the optimizer factory to extend the core EAGO optimizer for this special problem. The [`SubSolvers`](@ref) constructor is used to set the extension type (`t`), as well as the relaxed optimizer (`r`) and upper-bounding optimizer (`u`), if necessary. In this case, we will use the default solvers and only set the extension type.

```julia
factory = () -> Optimizer(SubSolvers(; t = QuasiConvex()))
```

## Construct the JuMP Model and Optimize

We now build the JuMP [[3](#References)] model representing this problem, solve it, and retrieve the solution.

```julia
opt = optimizer_with_attributes(factory, 
                                "absolute_tolerance" => 1E-8, 
                                "branch_variable" => branch_variable,
                                "iteration_limit" => 1000)
m = Model(opt)
@variable(m, ((i<6) ? 0.0 : -5.0) <= y[i=1:6] <= ((i<6) ? 5.0 : 0.0))
@constraint(m, sum(i*y[i] for i=1:5) - 5.0 == 0.0)
@constraint(m, sum(y[i]^2 for i=1:5) - 0.5*pi^2 <= 0.0)
@expression(m, expr1, 2.0*y[1]*y[2] + 4.0*y[1]*y[3] + 2.0*y[2]*y[3])
@constraint(m, -(0.5*y[1]^2 + 0.5*y[2]^2 + y[3]^2 + expr1) <= 0.0)
@NLexpression(m, expr2, log((5.0 + y[1])^2 + sum(y[i] for i=1:5)))
@NLconstraint(m, -y[1]^2 - 6.0*y[1]*y[2] - 2.0*y[2]^2 + cos(y[1]) + pi <= 0.0)
@NLconstraint(m, -expr2/(1.0 + sum(y[i]^2 for i=1:5)) - y[6] <= 0.0)
@objective(m, Min, y[6])

JuMP.optimize!(m)
```

## Retrieve Results

We then recover the solution values and the objective value using standard JuMP syntax.

```julia
solution = JuMP.value.(y[1:5])
global_obj_value = JuMP.value.(y[6])
print("Global solution at y*=$solution with a value of f*=$global_obj_value")
```

## References

1. C. Jansson, Quasiconvex relaxations based on interval arithmetic, Linear Algebra and its Applications, 324 (2001), pp. 27–53.
2. S. Boyd and L. Vandenberghe, Convex optimization, Cambridge University Press, 2004.
3. Iain Dunning and Joey Huchette and Miles Lubin. JuMP: A Modeling Language for Mathematical Optimization, *SIAM Review*, 59 (2017), pp. 295-320.
