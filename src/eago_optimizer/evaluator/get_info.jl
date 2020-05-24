# Copyright (c) 2018: Matthew Wilhelm & Matthew Stuber.
# This work is licensed under the Creative Commons Attribution-NonCommercial-
# ShareAlike 4.0 International License. To view a copy of this license, visit
# http://creativecommons.org/licenses/by-nc-sa/4.0/ or send a letter to Creative
# Commons, PO Box 1866, Mountain View, CA 94042, USA.
#############################################################################
# EAGO
# A development environment for robust and global optimization
# See https://github.com/PSORLab/EAGO.jl
#############################################################################
# src/eago_optimizer/evaluator/get_info.jl
# Access functions for information from evaluator.
#############################################################################

function MOI.eval_objective(d::Evaluator, x::Vector{Float64})
    forward_reverse_pass(d, x)
    val = 0.0
    if d.has_nlobj
        if d.objective.numvalued[1]
            val = d.objective.numberstorage[1]
        else
            val = d.objective.setstorage[1].cv
        end
    else
        error("No nonlinear objective.")
    end
    return val
end

"""
$(FUNCTIONNAME)

Retrieves the lower bound of the objective.
"""
function eval_objective_lo(d::Evaluator)
    val = 0.0
    if d.has_nlobj
        if d.objective.numvalued[1]
            val = d.objective.numberstorage[1]
        else
            val = d.objective.setstorage[1].Intv.lo
        end
    else
        error("No nonlinear objective.")
    end
    return val
end

function MOI.eval_constraint(d::Evaluator, g::Vector{Float64}, x::Vector{Float64})
    forward_reverse_pass(d, x)

    n = d.current_node

    print_info = false
    print_info && println("n = $n")
    #objvar = MC{11,NS}(x[1], Interval{Float64}(n.lower_variable_bounds[1], n.upper_variable_bounds[1]), 1)
    xMC = [MC{11,NS}(x[i], Interval{Float64}(n.lower_variable_bounds[i], n.upper_variable_bounds[i]), i) for i in 1:10]

    print_info && println("d.constraints[1].setstorage: $(length(d.constraints[1].nd))")
    for i = 1:length(d.constraints)
        if i === 1
            print_info && println(" ")
            print_info && println("i = $i, setstorage = $(d.constraints[i].setstorage[1])")
            #val = -xMC[1]*(1.12 + 0.13167*xMC[8] - 0.00667*(xMC[8])^2)+xMC[4] # interval and cc are different
            val9 = 0.0
            print_info && println("val9 = $(val9)")
            val8 = xMC[1]
            print_info && println("val8 = $(val8), type = $(typeof(val8))")
            val7 = xMC[8]
            print_info && println("val7 = $(val7)")
            val6 = 0.13167
            print_info && println("val6 = $(val6)")
            val5 = val6*val7
            print_info && println("val5 = $(val5)")
            val4 = 1.12
            print_info && println("val4 = $(val4)")
            val3 = val4 - val5
            print_info && println("val3 = $(val3), type = $(typeof(val3))")
            val2 = val8*val3
            print_info && println("val2 = $(val2), type = $(typeof(val2))")
            val1 = xMC[1]*(1.12 - 0.13167*xMC[8])
            print_info && println("val1 = $(val1)")

            #=
        elseif i === 2
            val = -0.001*xMC[4]*xMC[9]*xMC[6]/(98 - xMC[6]) + xMC[3]
        elseif i === 3
            #val = -(1.098*xMC[8] - 0.038*(xMC[8])^2) - 0.325*xMC[6] + xMC[7] - 57.425 # INTERVAL BOUNDS ARE DIFFERENT
            val = -1.098*xMC[8] + xMC[8]*xMC[8] - 0.325*xMC[6] + xMC[7] - 57.425
        elseif i === 4
            val = -(xMC[2] + xMC[5])/xMC[1] + xMC[8]
        elseif i === 5
            val = -0.063*xMC[4]*xMC[7] + 5.04*xMC[1] + 0.035*xMC[2] + 10*xMC[3] + 3.36*xMC[5] - objvar
            =#
        end
        #println("i = $i, val = $(val), val0 = $(val0)")
        print_info && println(" ")
        if d.constraints[i].numvalued[1]
            g[i] = d.constraints[i].numberstorage[1]
        else
            g[i] = d.constraints[i].setstorage[1].cv
        end
    end
    return
end

"""
$(FUNCTIONNAME)

Populates `g` with the concave relaxations of the constraints of `d` evaluated
at `y`.
"""
function eval_constraint_cc(d::Evaluator, g::Vector{Float64}, y::Vector{Float64})
    forward_reverse_pass(d,y)
    for i = 1:length(d.constraints)
        if d.constraints[i].numvalued[1]
            g[i] = d.constraints[i].numberstorage[1]
        else
            g[i] = d.constraints[i].setstorage[1].cc
        end
    end
    return
end

"""
$(FUNCTIONNAME)

Populates `g` with the lower bounds of the constraints of `d`.
"""
function eval_constraint_lo!(d::Evaluator, g::Vector{Float64})
    for i = 1:length(d.constraints)
        if d.constraints[i].numvalued[1]
            g[i] = d.constraints[i].numberstorage[1]
        else
            g[i] = d.constraints[i].setstorage[1].Intv.lo
        end
    end
    return
end

"""
$(FUNCTIONNAME)

Populates `g` with the upper bounds of the constraints of `d`.
"""
function eval_constraint_hi!(d::Evaluator, g::Vector{Float64})
    for i = 1:length(d.constraints)
        if d.constraints[i].numvalued[1]
            g[i] = d.constraints[i].numberstorage[1]
        else
            g[i] = d.constraints[i].setstorage[1].Intv.hi
        end
    end
    return
end


function MOI.eval_objective_gradient(d::Evaluator, df::Vector{Float64}, x::Vector{Float64})
    forward_reverse_pass(d,x)
    if d.has_nlobj
        if ~d.objective.numvalued[1]
            for j = 1:length(d.objective.setstorage[1].cv_grad)
                df[j] = d.objective.setstorage[1].cv_grad[j]
            end
        end
    else
        error("No nonlinear objective.")
    end
    return
end

function MOI.eval_constraint_jacobian(d::Evaluator,g,x)
    forward_reverse_pass(d,x)
    for i = 1:length(d.constraints)
        if ~d.constraints[i].numvalued[1]
            for j = 1:d.variable_number
                g[i,j] = d.constraints[i].setstorage[1].cv_grad[j]
            end
        else
            for j = 1:length(d.constraints[i].setstorage[1].cv_grad)
                g[i,j] = 0.0
            end
        end
    end
    return
end

"""
$(FUNCTIONNAME)

Populates `g` with the subgradients of the constraints of `d` evaluated at `y`.
"""
function eval_constraint_cc_grad(d::Evaluator, g, y)
        forward_reverse_pass(d,y)
        for i = 1:length(d.constraints)
            if ~d.constraints[i].numvalued[1]
                for j = 1:d.variable_number
                    g[i,j] = d.constraints[i].setstorage[1].cc_grad[j]
                end
            else
                for j = 1:d.variable_number
                    g[i,j] = 0.0
                end
            end
        end
    return
end

# looks good
function MOI.features_available(d::Evaluator)
    features = Symbol[]
    if !d.disable_1storder
        push!(features,:Grad)
        push!(features,:Jac)
    end
    return features
end

# looks good, doesn't do anything, EAGO builds the evaluator and attaches it to lower problems
function MOI.initialize(d::Evaluator, requested_features::Vector{Symbol})
end

# TO DO (CHECK GRADIENT DIMS)
function MOI.eval_constraint_jacobian_product(d::Evaluator, y, x, w)
    if !d.disable_1storder
        forward_reverse_pass(d,x)
        t = typeof(d.constraints[1].setstorage[1])
        y = zeros(t,length(d.constraints[1].setstorage[1].cv_grad),length(d.constraints))
        for i = 1:length(d.constraints)
            if ~d.constraints[i].numvalued[1]
                for j = 1:d.variable_number
                    y[i] += d.constraints[i].setstorage[1].cv_grad[j]*w[j]
                end
            end
        end
    else
        error("First order information unavailable.")
    end
    return
end

# TO DO
function MOI.eval_constraint_jacobian_transpose_product(d::Evaluator, y, x, w)
    if !d.disable_1storder
        forward_reverse_pass(d,x)
        y = zeros(Float64,length(d.constraints[1].setstorage[1].cv_grad),length(d.constraints))
        for i in 1:length(d.constraints)
            if ~d.constraints[i].numvalued[1]
                for j = 1:d.variable_number
                    y[i] += d.constraints[i].setstorage[1].cv_grad[j]*w[j]
                end
            end
        end
    else
        error("First order information unavailable.")
    end
    return
end

function MOI.jacobian_structure(d::Evaluator)
    jacobian_sparsity = Tuple{Int64,Int64}[]
    for row = 1:length(d.constraints)
        row_sparsity = d.constraints[row].grad_sparsity
        for idx in row_sparsity
            push!(jacobian_sparsity, (row, idx))
        end
    end
    return jacobian_sparsity
end

function grad_sparsity(d::Evaluator, j::Int64)
    sparsity = Int64[]
    if j == 1
        sparsity = d.objective.grad_sparsity
    else
        sparsity = d.constraints[j-1].grad_sparsity
    end
    return sparsity
end

#=
const SCALAR_ACCESS_FUNCTIONS = Union{cv, cc, lo, hi}
function eval_objective(::typeof{f}, d::NonlinearFunction{N}, x::Vector{Float64}) where {f <: SCALAR_ACCESS_FUNCTIONS, N}
    @assert d.has_nlobj "No nonlinear objective."
    forward_reverse_pass(d, x, 0)
    d.numvalued[1] ? d.numberstorage[1] : f(d.setstorage[1])::Float64
end
function eval_constraint(::typeof{f}, d::NonlinearFunction{N}, x::Vector{Float64}, i::Int) where {f <: SCALAR_ACCESS_FUNCTIONS, N}
    forward_reverse_pass(d, x, i)
    d.numvalued[1] ? d.numberstorage[1] : f(d.setstorage[1])
end

const VECTOR_ACCESS_FUNCTIONS = Union{cv_grad, cc_grad}
function eval_objective(::typeof{f}, d::NonlinearFunction{N}, df::Vector, x::Vector{Float64}) where {f <: VECTOR_ACCESS_FUNCTIONS, N}
    @assert d.has_nlobj "No nonlinear objective."
    forward_reverse_pass(d, x, 0)
    if d.numvalued[1]
        for i = 1:d.variable_number
            df[i] = f(d.setstorage[1])[i]
        end
    else
        fill!(df, 0.0)
    end
    nothing
end
function eval_objective(::typeof{f}, d::NonlinearFunction{N}, df::Vector, x::Vector{Float64}) where {f <: VECTOR_ACCESS_FUNCTIONS, N}
    forward_reverse_pass(d, x, 0)
    if d.numvalued[1]
        for i = 1:d.variable_number
            df[i] = f(d.setstorage[1])[i]
        end
    else
        fill!(df, 0.0)
    end
    nothing
end
grad_sparsity(d::NonlinearFunction) = d.grad_sparsity
=#
