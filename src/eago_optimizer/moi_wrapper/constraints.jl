MOI.supports_constraint(::Optimizer, ::Type{MOI.SingleVariable}, ::Type{MOI.LessThan{Float64}}) = true
MOI.supports_constraint(::Optimizer, ::Type{MOI.SingleVariable}, ::Type{MOI.GreaterThan{Float64}}) = true
MOI.supports_constraint(::Optimizer, ::Type{MOI.SingleVariable}, ::Type{MOI.EqualTo{Float64}}) = true
MOI.supports_constraint(::Optimizer, ::Type{MOI.SingleVariable}, ::Type{MOI.Interval{Float64}}) = true
MOI.supports_constraint(::Optimizer, ::Type{MOI.SingleVariable}, ::Type{MOI.ZeroOne}) = true

MOI.supports_constraint(::Optimizer, ::Type{MOI.ScalarAffineFunction{Float64}}, ::Type{MOI.LessThan{Float64}}) = true
MOI.supports_constraint(::Optimizer, ::Type{MOI.ScalarAffineFunction{Float64}}, ::Type{MOI.GreaterThan{Float64}}) = true
MOI.supports_constraint(::Optimizer, ::Type{MOI.ScalarAffineFunction{Float64}}, ::Type{MOI.EqualTo{Float64}}) = true
MOI.supports_constraint(::Optimizer, ::Type{MOI.ScalarAffineFunction{Float64}}, ::Type{MOI.Interval{Float64}}) = true

MOI.supports_constraint(::Optimizer, ::Type{MOI.ScalarQuadraticFunction{Float64}}, ::Type{MOI.LessThan{Float64}}) = true
MOI.supports_constraint(::Optimizer, ::Type{MOI.ScalarQuadraticFunction{Float64}}, ::Type{MOI.GreaterThan{Float64}}) = true
MOI.supports_constraint(::Optimizer, ::Type{MOI.ScalarQuadraticFunction{Float64}}, ::Type{MOI.EqualTo{Float64}}) = true
MOI.supports_constraint(::Optimizer, ::Type{MOI.ScalarQuadraticFunction{Float64}}, ::Type{MOI.Interval{Float64}}) = true

MOI.supports(::Optimizer, ::MOI.ObjectiveSense) = true

function MOI.add_variable(m::Optimizer)
    m.variable_number += 1
    m.nonlinear_variable[m.variable_number] = false
    m.fixed_variable[m.variable_number] = false
    push!(m.variable_info, VariableInfo())
    return MOI.VariableIndex(length(m.variable_info))
end
MOI.add_variables(m::Optimizer, n::Int) = [MOI.add_variable(m) for i in 1:n]


function MOI.add_constraint(m::Optimizer, v::MOI.SingleVariable, zo::MOI.ZeroOne)
    vi = v.variable
    check_inbounds(m, vi)
    has_upper_bound(m, vi) && error("Upper bound on variable $vi already exists.")
    has_lower_bound(m, vi) && error("Lower bound on variable $vi already exists.")
    is_fixed(m, vi) && error("Variable $vi is fixed. Cannot also set upper bound.")
    m.variable_info[vi.value].lower_bound = 0.0
    m.variable_info[vi.value].upper_bound = 1.0
    m.variable_info[vi.value].has_lower_bound = true
    m.variable_info[vi.value].has_upper_bound = true
    m.variable_info[vi.value].is_integer = true
    return MOI.ConstraintIndex{MOI.SingleVariable, MOI.ZeroOne}(vi.value)
end

function MOI.add_constraint(m::Optimizer, v::MOI.SingleVariable, lt::MOI.LessThan{Float64})
    vi = v.variable
    check_inbounds(m, vi)
    if isnan(lt.upper)
        error("Invalid upper bound value $(lt.upper).")
    end
    if has_upper_bound(m, vi)
        error("Upper bound on variable $vi already exists.")
    end
    if is_fixed(m, vi)
        error("Variable $vi is fixed. Cannot also set upper bound.")
    end
    m.variable_info[vi.value].upper_bound = lt.upper
    m.variable_info[vi.value].has_upper_bound = true
    return MOI.ConstraintIndex{MOI.SingleVariable, MOI.LessThan{Float64}}(vi.value)
end

function MOI.add_constraint(m::Optimizer, v::MOI.SingleVariable, gt::MOI.GreaterThan{Float64})
    vi = v.variable
    check_inbounds(m, vi)
    if isnan(gt.lower)
        error("Invalid lower bound value $(gt.lower).")
    end
    if has_lower_bound(m, vi)
        error("Lower bound on variable $vi already exists.")
    end
    if is_fixed(m, vi)
        error("Variable $vi is fixed. Cannot also set lower bound.")
    end
    m.variable_info[vi.value].lower_bound = gt.lower
    m.variable_info[vi.value].has_lower_bound = true
    return MOI.ConstraintIndex{MOI.SingleVariable, MOI.GreaterThan{Float64}}(vi.value)
end

function MOI.add_constraint(m::Optimizer, v::MOI.SingleVariable, eq::MOI.EqualTo{Float64})
    vi = v.variable
    check_inbounds(m, vi)
    if isnan(eq.value)
        error("Invalid fixed value $(gt.lower).")
    end
    if has_lower_bound(m, vi)
        error("Variable $vi has a lower bound. Cannot be fixed.")
    end
    if has_upper_bound(m, vi)
        error("Variable $vi has an upper bound. Cannot be fixed.")
    end
    if is_fixed(m, vi)
        error("Variable $vi is already fixed.")
    end
    m.variable_info[vi.value].lower_bound = eq.value
    m.variable_info[vi.value].upper_bound = eq.value
    m.variable_info[vi.value].has_lower_bound = true
    m.variable_info[vi.value].has_upper_bound = true
    m.variable_info[vi.value].is_fixed = true
    return MOI.ConstraintIndex{MOI.SingleVariable, MOI.EqualTo{Float64}}(vi.value)
end

function MOI.add_constraint(m::Optimizer, v::MOI.SingleVariable, eq::MOI.Interval{Float64})
    vi = v.variable
    check_inbounds(m, vi)
    if isnan(eq.lower)
        error("Invalid fixed value $(gt.lower).")
    end
    if isnan(eq.upper)
        error("Invalid fixed value $(gt.upper).")
    end
    if has_lower_bound(m, vi)
        error("Lower bound on variable $vi already exists. Cannot also set interval bounds.")
    end
    if has_upper_bound(m, vi)
        error("Upper bound on variable $vi already exists. Cannot also set interval bounds.")
    end
    if is_fixed(m, vi)
        error("Variable $vi is fixed. Cannot also set interval bounds.")
    end
    m.variable_info[vi.value].lower_bound = eq.lower
    m.variable_info[vi.value].upper_bound = eq.upper
    m.variable_info[vi.value].has_lower_bound = true
    m.variable_info[vi.value].has_upper_bound = true
    return MOI.ConstraintIndex{MOI.SingleVariable, MOI.Interval{Float64}}(vi.value)
end

macro define_addconstraint_linear(function_type, set_type, array_name)
    quote
        function MOI.add_constraint(m::Optimizer, func::$function_type, set::$set_type)
            check_inbounds(m, func)
            push!(m.$(array_name), (func, set, length(func.terms)))
            indx = MOI.ConstraintIndex{$function_type, $set_type}(length(m.$(array_name)))
            m.constraint_convexity[indx] = true
            return indx
        end
    end
end

macro define_addconstraint_quadratic(function_type, set_type, array_name)
    quote
        function MOI.add_constraint(m::Optimizer, func::$function_type, set::$set_type)
            check_inbounds(m, func)
            for i in func.affine_terms m.nonlinear_variable[i.variable_index.value] = true end
            for i in func.quadratic_terms
                m.nonlinear_variable[i.variable_index_1.value] = true
                m.nonlinear_variable[i.variable_index_2.value] = true
            end
            push!(m.$(array_name), (func, set, length(m.$(array_name))+1))
            indx = MOI.ConstraintIndex{$function_type, $set_type}(length(m.$(array_name)))
            m.constraint_convexity[indx] = false
            return indx
        end
    end
end

@define_addconstraint_linear MOI.ScalarAffineFunction{Float64} MOI.LessThan{Float64} linear_leq_constraints
@define_addconstraint_linear MOI.ScalarAffineFunction{Float64} MOI.GreaterThan{Float64} linear_geq_constraints
@define_addconstraint_linear MOI.ScalarAffineFunction{Float64} MOI.EqualTo{Float64} linear_eq_constraints
@define_addconstraint_linear MOI.ScalarAffineFunction{Float64} MOI.Interval{Float64} linear_interval_constraints

@define_addconstraint_quadratic MOI.ScalarQuadraticFunction{Float64} MOI.LessThan{Float64} quadratic_leq_constraints
@define_addconstraint_quadratic MOI.ScalarQuadraticFunction{Float64} MOI.GreaterThan{Float64} quadratic_geq_constraints
@define_addconstraint_quadratic MOI.ScalarQuadraticFunction{Float64} MOI.EqualTo{Float64} quadratic_eq_constraints
@define_addconstraint_quadratic MOI.ScalarQuadraticFunction{Float64} MOI.Interval{Float64} quadratic_interval_constraints

function MOI.set(m::Optimizer, ::MOI.NLPBlock, nlp_data::MOI.NLPBlockData)
    m.nlp_data = nlp_data
    return
end
