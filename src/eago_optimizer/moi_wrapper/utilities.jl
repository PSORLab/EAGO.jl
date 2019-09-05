function MOI.empty!(m::Optimizer)
    m = Optimizer()
end

function MOI.is_empty(m::Optimizer)

    vb = Bool[false for i=1:3]

    vb[1] = isempty(m.variable_info)
    vb[2] = m.optimization_sense == MOI.FEASIBILITY_SENSE
    vb[3] = m.termination_status_code == MOI.OPTIMIZE_NOT_CALLED

    bool_out = (sum(vb) > 0)

    return bool_out
end

function check_inbounds(m::Optimizer, vi::MOI.VariableIndex)
    num_variables = length(m.variable_info)
    if !(1 <= vi.value <= num_variables)
        error("Invalid variable index $vi. ($num_variables variables in the model.)")
    end
end

check_inbounds(m::Optimizer, var::MOI.SingleVariable) = check_inbounds(m, var.variable)

function check_inbounds(m::Optimizer, aff::MOI.ScalarAffineFunction)
    for term in aff.terms
        check_inbounds(m, term.variable_index)
    end
end

function check_inbounds(m::Optimizer, quad::MOI.ScalarQuadraticFunction)
    for term in quad.affine_terms
        check_inbounds(m, term.variable_index)
    end
    for term in quad.quadratic_terms
        check_inbounds(m, term.variable_index_1)
        check_inbounds(m, term.variable_index_2)
    end
end

function has_upper_bound(m::Optimizer, vi::MOI.VariableIndex)
    return m.variable_info[vi.value].has_upper_bound
end

function has_lower_bound(m::Optimizer, vi::MOI.VariableIndex)
    return m.variable_info[vi.value].has_lower_bound
end

function is_fixed(m::Optimizer, vi::MOI.VariableIndex)
    return m.variable_info[vi.value].is_fixed
end

function is_integer_feasible(m::Optimizer)
    int_feas = true
    for var in m.integer_variables
        (0.0 < m.current_lower_info.solution[var] < 1.0) && (int_feas = false; break)
    end
    return int_feas
end

is_integer_variable(m::Optimizer, i::Int) = in(i, m.integer_variables)

function ReverseDict(dict)
    Dict(value => key for (key, value) in dict)
end

function MOI.copy_to(model::Optimizer, src::MOI.ModelLike; copy_names = false)
    return MOI.Utilities.default_copy_to(model, src, copy_names)
end

"""
    create_initial_node!

Creates an initial node with initial box constraints and adds it to the stack
"""
function create_initial_node!(m::Optimizer)
    m.stack[1] = NodeBB()
    m.stack[1].lower_variable_bounds = lower_bound.(m.variable_info)
    m.stack[1].upper_variable_bounds = upper_bound.(m.variable_info)
    m.current_node_count = 1
    m.maximum_node_id += 1
end

function label_nonlinear_variables!(m::Optimizer, x)
    # scans subexpressions, objective, and constraints for nonlinear terms
    if ~isa(x, EmptyNLPEvaluator)
        if x.has_nlobj
            if (x.objective.linearity != JuMP._Derivatives.LINEAR) ||
               (x.objective.linearity != JuMP._Derivatives.CONSTANT)
                for i in 1:length(x.objective.nd)
                    nd = x.objective.nd[i]
                    if (nd.nodetype == JuMP._Derivatives.VARIABLE)
                        m.bisection_variable[nd.index] = true
                    end
                end
            end
        end
        for i in 1:length(x.constraints)
            if (x.constraints[i].linearity != JuMP._Derivatives.LINEAR) ||
               (x.constraints[i].linearity != JuMP._Derivatives.CONSTANT)
                for j in 1:length(x.constraints[i].nd)
                    nd = x.constraints[i].nd[j]
                    bool1 = (nd.nodetype == JuMP._Derivatives.VARIABLE)
                    if (nd.nodetype == JuMP._Derivatives.VARIABLE)
                        m.bisection_variable[nd.index] = true
                    end
                end
            end
        end
        for i in 1:length(x.subexpressions)
            if (x.subexpressions[i].linearity != JuMP._Derivatives.LINEAR) ||
               (x.subexpressions[i].linearity != JuMP._Derivatives.CONSTANT)
                for j in 1:length(x.subexpressions[i].nd)
                    nd = x.subexpressions[i].nd[j]
                    bool1 = (nd.nodetype == JuMP._Derivatives.VARIABLE)
                    if (nd.nodetype == JuMP._Derivatives.VARIABLE)
                        m.bisection_variable[nd.index] = true
                    end
                end
            end
        end
    end
end

function get_objective_bounds(m::Optimizer)
end

function add_optimization_variable()
end
