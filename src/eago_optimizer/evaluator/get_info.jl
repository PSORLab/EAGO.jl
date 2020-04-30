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

#=
function MOI.eval_objective(d::Evaluator, x::Vector{Float64})
    forward_reverse_pass(d,x)
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
    forward_reverse_pass(d,x)
    for i = 1:length(d.constraints)
        #println("i = $i, setstorage = $(d.constraints[i].setstorage[1])")
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
            for j in 1:length(d.objective.setstorage[1].cv_grad)
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
            for j in 1:d.variable_number
                g[i,j] = d.constraints[i].setstorage[1].cv_grad[j]
            end
        else
            for j in 1:length(d.constraints[i].setstorage[1].cv_grad)
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
                for j in 1:d.variable_number
                    g[i,j] = d.constraints[i].setstorage[1].cc_grad[j]
                end
            else
                for j in 1:d.variable_number
                    g[i,j] = 0.0
                end
            end
        end
    return
end
=#
# looks good
#=
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

function grad_sparsity(d::Evaluator, j::Int64)
    sparsity = Int64[]
    if j == 1
        sparsity = d.objective.grad_sparsity
    else
        sparsity = d.constraints[j-1].grad_sparsity
    end
    return sparsity
end
=#

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
