
###
### AFFINE FUNCTIONS
###

function lower_interval_bound(m::Optimizer, f::AffineFunctionIneq, y::NodeBB)

    terms = f.terms
    lo_bnds = y.lower_variable_bounds
    up_bnds = y.upper_variable_bounds

    lower_interval_bound = f.constant
    for i = 1:f.len
        coeff, indx = @inbounds terms[i]

        if m._branch_variables[vi]
            mapped_vi = @inbounds sol_branch_map[indx]
            xL = @inbounds lo_bnds[mapped_vi]
            xU = @inbounds up_bnds[mapped_vi]
        else
            xL = @inbounds m._working_problem._variable_info[indx].lower_bound
            xU = @inbounds m._working_problem._variable_info[indx].upper_bound
        end

        lower_interval_bound += (coeff > 0.0) ? coeff*xL : coeff*xU
    end

    return lower_interval_bound
end

function interval_bound(m::Optimizer, f::AffineFunctionEq, y::NodeBB)
    terms = f.terms
    lo_bnds = y.lower_variable_bounds
    up_bnds = y.upper_variable_bounds

    lower_interval_bound = f.constant
    upper_interval_bound = f.constant
    for i = 1:f.len

        if m._branch_variables[vi]
            mapped_vi = @inbounds sol_branch_map[indx]
            xL = @inbounds lo_bnds[mapped_vi]
            xU = @inbounds up_bnds[mapped_vi]
        else
            xL = @inbounds m._working_problem._variable_info[indx].lower_bound
            xU = @inbounds m._working_problem._variable_info[indx].upper_bound
        end

        if coeff > 0.0
            lower_interval_bound += coeff*xL
            upper_interval_bound += coeff*xU
        else
            lower_interval_bound += coeff*xU
            upper_interval_bound += coeff*xL
        end
    end

    return lower_interval_bound, upper_interval_bound
end

###
### QUADRATIC FUNCTIONS
###

function lower_interval_bound(m::Optimizer, f::BufferedQuadraticIneq, n::NodeBB)

    sol_branch_map = m._sol_to_branch_map
    lo_bnds = n.lower_variable_bounds
    up_bnds = n.upper_variable_bounds
    lower_interval_bound = Interval{Float64}(f.func.constant)

    for aff_term in f.func.affine_terms

        coeff = aff_term.coefficient

        vi = aff_term.variable_index.value
        if m._branch_variables[vi]
            mapped_vi = @inbounds sol_branch_map[vi]
            xL = @inbounds lo_bnds[mapped_vi]
            xU = @inbounds up_bnds[mapped_vi]
        else
            xL = @inbounds m._working_problem._variable_info[vi].lower_bound
            xU = @inbounds m._working_problem._variable_info[vi].upper_bound
        end

        lower_interval_bound += coeff > 0.0 ? coeff*xL : coeff*xU
    end

    # assumes branching on all quadratic terms, otherwise we'd need to distinguish
    # between look ups to the node bounds and lookups to the variable info in the
    # working problem.
    for quad_term in f.func.quadratic_terms

        coeff = quad_term.coefficient

        vi1 = quad_term.variable_index_1.value
        vi2 = quad_term.variable_index_2.value
        mapped_vi1 = sol_branch_map[vi1]

        xL = @inbounds lo_bnds[mapped_vi1]
        xU = @inbounds up_bnds[mapped_vi1]

        if vi1 === vi2
            if coeff > 0.0
                lower_interval_bound += (0.0 < xL) ? 0.5*coeff*xL*xL : ((xU <= 0.0) ? 0.5*coeff*xU*xU : 0.0)
            else
                lower_interval_bound += (xL < xU) ? 0.5*coeff*xU*xU : 0.5*coeff*xL*xL
            end
        else
            mapped_vi2 = sol_branch_map[vi2]

            il2b = @inbounds lo_bnds[mapped_vi2]
            iu2b = @inbounds up_bnds[mapped_vi2]
            lower_interval_bound += coeff*Interval{Float64}(xL, xU)*Interval{Float64}(il2b, iu2b)

        end
    end

    return lower_interval_bound.lo
end

function interval_bound(m::Optimizer, f::BufferedQuadraticEq, n::NodeBB)

    lo_bnds = n.lower_variable_bounds
    up_bnds = n.upper_variable_bounds
    val_intv = Interval(f.func.constant)

    for aff_term in f.func.affine_terms

        coeff = aff_term.coefficient

        vi = aff_term.variable_index.value
        if m._branch_variables[vi]
            mapped_vi = @inbounds sol_branch_map[vi]
            xL = @inbounds lo_bnds[mapped_vi]
            xU = @inbounds up_bnds[mapped_vi]
        else
            xL = @inbounds m._working_problem._variable_info[vi].lower_bound
            xU = @inbounds m._working_problem._variable_info[vi].upper_bound
        end

        val_intv += coeff*Interval(xL, xU)
    end

    # assumes branching on all quadratic terms, otherwise we'd need to distinguish
    # between look ups to the node bounds and lookups to the variable info in the
    # working problem.
    for quad_term in f.func.quadratic_terms

        coeff = quad_term.coefficient

        vi1 = quad_term.variable_index_1.value
        vi2 = quad_term.variable_index_2.value

        mapped_vi1 = @inbounds sol_branch_map[vi1]

        xL = @inbounds lo_bnds[mapped_vi1]
        xU = @inbounds up_bnds[mapped_vi1]

        if vi1 === vi2
            val_intv += 0.5*coeff*pow(Interval(xL, xU), 2)
        else
            mapped_vi2 = @inbounds sol_branch_map[vi2]
            @inbounds il2b = lo_bnds[mapped_vi2]
            @inbounds iu2b = up_bnds[mapped_vi2]
            val_intv += coeff*Interval(xL, xU)*Interval(il2b, iu2b)
        end
    end

    return val_intv.lo, val_intv.hi
end

###
### NONLINEAR FUNCTIONS
###

function lower_interval_bound(m::Optimizer, d::BufferedNonlinearFunction{V}, n::NodeBB) where V
    #println("d.has_value: $(d.has_value)")
    if !d.has_value
        forward_pass!(m._working_problem._relaxed_evaluator, d)
    end

    expr = d.expr
    if expr.isnumber[1]
        lower_value = expr.numberstorage[1]
    else
        lower_value = expr.setstorage[1].Intv.lo
    end

    return lower_value
end

function interval_bound(m::Optimizer, d::BufferedNonlinearFunction{V}, n::NodeBB) where V
    if !d.has_value
        forward_pass!(d.evaluator, d)
    end

    expr = d.expr
    if expr.isnumber[1]
        interval_value = Interval(expr.numberstorage[1])
    else
        interval_value = expr.setstorage[1].Intv
    end

    return interval_value
end
