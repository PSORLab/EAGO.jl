#Runs the Poor Man's LP to refine variables present nonlinear using linear constraints
function poor_man_lp(m::Optimizer,n::NodeBB)

    # Runs Poor Man LP on constraints of form ax >= b
    for (func,constr,ind) in m.linear_geq_constraints
        TempValue = (constr.lower - func.constant)
        for term in func.terms
            ti = m.variable_index_to_storage[term.variable_index.value]
            TempValue += -max(term.coefficient*n.upper_variable_bounds[ti], term.coefficient*n.lower_variable_bounds[ti])
        end
        for term in func.terms
            vi = m.variable_index_to_storage[term.variable_index.value]
            TermValue = -max(term.coefficient*n.upper_variable_bounds[vi], term.coefficient*n.lower_variable_bounds[vi])
            CutValue = (TempValue - TermValue)/term.coefficient
            if (term.coefficient > 0.0 )
                if (n.lower_variable_bounds[vi] < CutValue)
                    (CutValue > n.upper_variable_bounds[vi]) && (return false)
                    n.lower_variable_bounds[vi] = CutValue
                end
            else
                if (n.upper_variable_bounds[vi] > CutValue)
                    (CutValue < n.lower_variable_bounds[vi]) && (return false)
                    n.upper_variable_bounds[vi] = CutValue
                end
            end
        end
    end
        # Runs Poor Man LP on constraints of form ax <= b
    for (func,constr,ind) in m.linear_leq_constraints
        TempValue = (constr.upper - func.constant)
        for term in func.terms
            ti = m.variable_index_to_storage[term.variable_index.value]
            TempValue += -min(term.coefficient*n.upper_variable_bounds[ti], term.coefficient*n.lower_variable_bounds[ti])
        end
        for term in func.terms
            vi = m.variable_index_to_storage[term.variable_index.value]
            TermValue = -min(term.coefficient*n.upper_variable_bounds[vi], term.coefficient*n.lower_variable_bounds[vi])
            CutValue = (TempValue - TermValue)/term.coefficient
            if (term.coefficient < 0.0 )
                if (n.lower_variable_bounds[vi] < CutValue)
                    (CutValue > n.upper_variable_bounds[vi]) && (return false)
                    n.lower_variable_bounds[vi] = CutValue
                end
            else
                if (n.upper_variable_bounds[vi] > CutValue)
                    (CutValue < n.lower_variable_bounds[vi]) && (return false)
                    n.upper_variable_bounds[vi] = CutValue
                end
            end
        end
    end

    for (func,constr,ind) in m.linear_eq_constraints
        TempValue = (constr.value - func.constant)
        for term in func.terms
            ti = m.variable_index_to_storage[term.variable_index.value]
            TempValue += -max(term.coefficient*n.upper_variable_bounds[ti], term.coefficient*n.lower_variable_bounds[ti])
        end
        for term in func.terms
            vi = m.variable_index_to_storage[term.variable_index.value]
            TermValue = -max(term.coefficient*n.upper_variable_bounds[vi], term.coefficient*n.lower_variable_bounds[vi])
            CutValue = (TempValue - TermValue)/term.coefficient
            if (term.coefficient > 0.0 )
                if (n.lower_variable_bounds[vi] < CutValue)
                    (CutValue > n.upper_variable_bounds[vi]) && (return false)
                    n.lower_variable_bounds[vi] = CutValue
                end
            else
                if (n.upper_variable_bounds[vi] > CutValue)
                    (CutValue < n.lower_variable_bounds[vi]) && (return false)
                    n.upper_variable_bounds[vi] = CutValue
                end
            end
        end
        TempValue = (constr.value - func.constant)
        for term in func.terms
            ti = m.variable_index_to_storage[term.variable_index.value]
            TempValue += -min(term.coefficient*n.upper_variable_bounds[ti], term.coefficient*n.lower_variable_bounds[ti])
        end
        for term in func.terms
            vi = m.variable_index_to_storage[term.variable_index.value]
            TermValue = -min(term.coefficient*n.upper_variable_bounds[vi], term.coefficient*n.lower_variable_bounds[vi])
            CutValue = (TempValue - TermValue)/term.coefficient
            if (term.coefficient < 0.0 )
                if (n.lower_variable_bounds[vi] < CutValue)
                    (CutValue > n.upper_variable_bounds[vi]) && (return false)
                    n.lower_variable_bounds[vi] = CutValue
                end
            else
                if (n.upper_variable_bounds[vi] > CutValue)
                    (CutValue < n.lower_variable_bounds[vi]) && (return false)
                    n.upper_variable_bounds[vi] = CutValue
                end
            end
        end
    end

    return true
end
