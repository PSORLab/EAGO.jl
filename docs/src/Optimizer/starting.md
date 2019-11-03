# Solving simple example: an artificial neural network with EAGO

[Matthew Wilhelm](https://psor.uconn.edu/person/matthew-wilhelm/)  
Department of Chemical and Biomolecular Engineering, University of Connecticut

In [1,2], a surrogate ANN model of bioreactor productivity was constructed by fitting results from computationally expensive CFD simulations. The author then optimized this surrogate model to obtain ideal processing conditions. This optimization problem is given by:

$
\begin{align}
    \label{prob:ANN1}
    \max_{\mathbf x\in X} B_2 + \sum_{i=1}^{3}\frac{2D_{i}}{1+\exp(-2y_{i})} \qquad \text{where} \qquad y_{i} =  B_i + \sum_{i=1}^{3}\sum_{j=1}^{8} W_{ij}x_{i}\\
\end{align}
$

We'll repeat this exercise using JuMP and EAGO.

### Input parameters

In the first block, we input parameters values supplied in the paper for $W$, $B$, $D$, and $X$ into Julia as simple array objects.

```julia
using JuMP, EAGO

# Box constraints for input variables
xLBD = [0.623   0.093   0.259   6.56   1114   0.013   0.127   0.004]
xUBD = [5.89    0.5     1.0     90     25000  0.149   0.889   0.049]

# Weights associated with the hidden layer
W = [ 0.54  -1.97  0.09  -2.14  1.01  -0.58  0.45  0.26;
     -0.81  -0.74  0.63  -1.60 -0.56  -1.05  1.23  0.93;
     -0.11  -0.38 -1.19   0.43  1.21   2.78 -0.06  0.40]

# Weights associated with the output layer
D = [-0.91 0.11 0.52]

# Bias associated with the hidden layer
B1 = [-2.698 0.012 2.926]

# Bias associated with the output layer
B2 = -0.46
```

### Construct the JuMP model and optimize

We now formulate the problem using standard JuMP[3] syntax and optimize it. Note that we are forming an NLexpression object to handle the summation term to keep the code visually simple but this could be placed directly in the JuMP expressions instead.

```julia
# Model construction
model = Model(with_optimizer(EAGO.Optimizer, absolute_tolerance = 0.001))
@variable(model, xLBD[i] <= x[i=1:8] <= xUBD[i])
@NLexpression(model, prop[i=1:3], B1[i] + sum(W[i,j]*x[i] for j in 1:8))
@NLobjective(model, Max, B2 + sum(D[i]*(2/(1+exp(-2*prop[i]))) for i=1:3))

# Solves the model
optimize!(model)
```

### Retrieve results

We then recover the objective value, the solution value, and termination status codes using standard JuMP syntax.

```julia
# Access calculated values
fval = JuMP.objective_value(model)
xsol = JuMP.value.(x)
status_term = JuMP.termination_status(model)
status_prim = JuMP.primal_status(model)

println("EAGO terminated with a status of $status_term and a result code of $status_prim")
println("The optimal value is: $fval, the solution found is $xsol.")
```

### Reference:
1. J. D. Smith, A. A. Neto, S. Cremaschi, and D. W. Crunkleton, CFD-based optimization of a flooded bed algae bioreactor, *Industrial & Engineering Chemistry Research*, 52 (2012), pp. 7181–7188
2. A. M. Schweidtmann and A. Mitsos. Global Deterministic Optimization with Artificial Neural Networks Embedded [https://arxiv.org/pdf/1801.07114.pdf](https://arxiv.org/pdf/1801.07114.pdf)
3. Iain Dunning and Joey Huchette and Miles Lubin. JuMP: A Modeling Language for Mathematical Optimization, *SIAM Review*, 59 (2017), pp. 295-320.
