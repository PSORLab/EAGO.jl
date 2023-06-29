# Advanced-Use Example 2

## An ``\alpha``BB Example for a QCQP.

In this example, we will demonstrate the use of a user-defined lower-bounding problem that uses ``\alpha``BB convex relaxations. In this example, we wish to solve the (nonconvex) QCQP:

```math
\begin{aligned}
&\min_{\mathbf x\in\mathbb R^2}\frac{1}{2}\mathbf x^{\rm T}\mathbf Q_f\mathbf x+\mathbf c_f^{\rm T}\mathbf x\\
{\rm s.t.}\;\;&g_1(\mathbf x)=\frac{1}{2}\mathbf x^{\rm T}\mathbf Q_{g_1}\mathbf x+\mathbf c_{g_1}^{\rm T}\mathbf x\le 0\\
&g_2(\mathbf x)=\frac{1}{2}\mathbf x^{\rm T}\mathbf Q_{g_2}\mathbf x+\mathbf c_{g_2}^{\rm T}\mathbf x\le 0\\
\end{aligned}
```

with ``\mathbf Q_i\in\mathbb R^{2\times 2}`` not positive semidefinite for any ``i``.

```julia
using JuMP, EAGO, Ipopt, LinearAlgebra
```

## Custom Function Definitions

For convenience, we'll define the following function that returns all the problem data ``\mathbf Q_i`` and ``\mathbf c_i``.

```julia
function QCQP_setup()

    Qf = [3. 3/2.;3/2. -5.]
    cf = [3.;2.]

    Qg1 = [-2. 5.;5. -2.]
    cg1 = [1.;3.]

    Qg2 = [-6. 3.;3. 2.]
    cg2 = [2.;1.]
    
    return Qf,cf,Qg1,cg1,Qg2,cg2
end
```

The next function we'll define will take as input data for a particular quadratic function and the interval bounds on the decision variables, and construct an ``\alpha``BB convex relaxation of that function. Since we're solving a QCQP, we'll use the ``\verb|eigvals|`` function to directly compute the eigenvalues of the input ``\mathbf Q_i`` matrix.

```julia
function αBB_relax(Q::Matrix{T},c::Vector{T},xL::Vector{T},xU::Vector{T},x::Real...) where {T<:Float64}
    α=max(0,-minimum(eigvals(Q))/2)
    y = [x[1];x[2]]
    cv = 1/2*y'*Q*y+c'*y+α*(xL-y)'*(xU-y)
    return cv
end
```

## Redefine the Lower-Bounding Problem

The following code first defines our EAGO extension (custom version) struct and then it redefines the lower-bounding problem as our own version. That is, when we call this customized version of EAGO to solve the problem, it'll deploy this version of the lower-bounding problem instead of the default version.  

```julia
import EAGO: Optimizer, GlobalOptimizer

struct αBB_Convex <: EAGO.ExtensionType end
import EAGO: lower_problem!
function EAGO.lower_problem!(t::αBB_Convex, opt::GlobalOptimizer)
    # Get active node
    n = opt._current_node
    # Get bounds on active node for calculating relaxations
    xL = n.lower_variable_bounds[1:2]
    xU = n.upper_variable_bounds[1:2]
    # Get the problem data
    Qf,cf,Qg1,cg1,Qg2,cg2=QCQP_setup()

    # Define the JuMP model and declare the solver
    mL = JuMP.Model(JuMP.optimizer_with_attributes(Ipopt.Optimizer,
                                "tol" => 1.0e-6,
                                "print_level" => 0))
    @variable(mL,xL[i]<=x[i=1:2]<=xU[i])
    
    # Define the function closures for the user-defined relaxations
    fcv(x...) = αBB_relax(Qf,cf,xL,xU,x...)
    g1cv(x...) = αBB_relax(Qg1,cg1,xL,xU,x...)
    g2cv(x...) = αBB_relax(Qg2,cg2,xL,xU,x...)

    # Register the user-defined functions
    # Note: if the gradients and Hessians are directly available, they could
    # be passed as arguments to the register function to speed things up.
    JuMP.register(mL,:fcv,2,fcv,autodiff=true)
    JuMP.register(mL,:g1cv,2,g1cv,autodiff=true)
    JuMP.register(mL,:g2cv,2,g2cv,autodiff=true)

    # Declare the objective function and constraints
    @NLobjective(mL,Min,fcv(x[1],x[2]))
    @NLconstraint(mL,g1cv(x[1],x[2])<=0.)
    @NLconstraint(mL,g2cv(x[1],x[2])<=0.)
    
    # Solve the relaxed problem
    JuMP.optimize!(mL)
    
    # Get primal status, termination status, determine if a global solution was obtained
    tstatus = MOI.get(mL, MOI.TerminationStatus())
    pstatus = MOI.get(mL, MOI.PrimalStatus())

    solution = JuMP.value.(x)
    # Interpret status codes for branch-and-bound
    if EAGO.local_problem_status(tstatus, pstatus) == EAGO.LRS_FEASIBLE
        opt._lower_objective_value = JuMP.objective_value(mL) 
        opt._lower_solution[1:length(solution)] = solution
        opt._lower_feasibility = true
        opt._cut_add_flag = false
    else
        opt._lower_feasibility = false
        opt._lower_objective_value = -Inf
        opt._cut_add_flag = false
    end
    return
end
```

!!! note 

    Caution: By default, EAGO solves the epigraph reformulation of your original problem, which increases the original problem dimensionality by +1 with the introduction of an auxiliary variable. When defining custom routines (such as the lower-bounding problem here) that are intended to work nicely with default EAGO routines (such as preprocessing), the user must account for the *new* dimensionality of the problem. In the code above, we wish to access the information of the specific B&B node and define an optimization problem based on that information. However, in this example, the node has information for 3 variables (the original 2 plus 1 for the auxiliary variable appended to the original variable vector) as ``(x_1,x_2,\eta)``. The lower-bounding problem was defined to optimize the relaxed problem with respect to the original 2 decision variables. When storing the results of this subproblem to the current B&B node, it is important to take care to store the information at the appropriate indices and not inadvertently redefine the problem dimensionality (i.e., by simply storing the optimization solution as the `lower_solution` of the current node). For problems that are defined to only branch on a subset of the original variables, the optimizer has a member `_sol_to_branch_map` that carries the mapping between the indices of the original variables to those of the variables being branched on. See the custom_quasiconvex example to see how this is done.


## Turn Off Preprocessing and Postprocessing Routines

(Optional) Turn off (short circuit) preprocessing routines if you don't want to use them as defined in EAGO. 

```julia
import EAGO: preprocess!
function EAGO.preprocess!(t::αBB_Convex, x::GlobalOptimizer)
    x._preprocess_feasibility = true
    return
end
```

(Optional) Turn off (short circuit) postprocessing routines if you don't want to use them as defined in EAGO. 

```julia
import EAGO: postprocess!
function EAGO.postprocess!(t::αBB_Convex, x::GlobalOptimizer)
    x._postprocess_feasibility = true
    return
end
```

## Construct the JuMP Model and Optimize

Now, we'll tell EAGO to use our custom/extended solver, set up the main JuMP model, and solve it with our custom solver. 

```julia
factory = () -> EAGO.Optimizer(SubSolvers(; t=αBB_Convex() ))
m = JuMP.Model(optimizer_with_attributes(factory,
                                "relative_tolerance" => 1e-3,
                                "verbosity" => 1,
                                "output_iterations" => 1, 
                                "branch_variable" => Bool[true; true],
                                ))
Qf,cf,Qg1,cg1,Qg2,cg2=QCQP_setup()  # Get QCQP data
xL = [-3.;-5.]                      # Lower bounds on x
xU = [1.; 2]                        # Upper bounds on x
@variable(m, xL[i] <= x[i=1:2] <= xU[i])

# Define objective and constraints
@objective(m,Min,1/2*x'*Qf*x+cf'*x)
@constraint(m,1/2*x'*Qg1*x+cg1'*x<=0.)
@constraint(m,1/2*x'*Qg2*x+cg2'*x<=0.)

# Solve the problem
@time optimize!(m)
```

## Retrieve Results

We then recover the objective value, the solution value, and termination status codes using standard JuMP syntax.

```julia
println("x1* = ", JuMP.value(x[1]), " x2* = ",
         JuMP.value(x[2])," f* = ",JuMP.objective_value(m))
TermStatus = JuMP.termination_status(m)
PrimStatus = JuMP.primal_status(m)
println("Algorithm terminated with a status of $TermStatus and result code of $PrimStatus")
```