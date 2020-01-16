"""
$(FUNCTIONNAME)

Tightens the bounds of the `_current_node` using the current global upper bound
and the duality information obtained from the relaxation.
"""
function variable_dbbt!(x::NodeBB, mult_lo::Vector{Float64}, mult_hi::Vector{Float64},
                        LBD::Float64, UBD::Float64, nx::Int64)

    sempty = isempty(x)
    oldx = deepcopy(x)
    cut = 0.0
    vb = 0.0
    lvbs = x.lower_variable_bounds
    uvbs = x.upper_variable_bounds
    if LBD <= UBD
        for i = 1:nx
            @inbounds ml = mult_lo[i]
            @inbounds mh = mult_hi[i]
            if ml > 0.0
                @inbounds cut = lvbs[i] + (UBD - LBD)/ml
                @inbounds vb = uvbs[i]
                if cut < vb
                    @inbounds uvbs[i] = cut
                end
            elseif mh > 0.0
                 @inbounds cut = uvbs[i] - (UBD - LBD)/mh
                 @inbounds vb = lvbs[i]
                 if cut > vb
                     @inbounds lvbs[i] = cut
                 end
            end
         end
    end
    return
end

"""
$(FUNCTIONNAME)

Excludes OBBT on variable indices that are tight for the solution of the relaxation.
"""
function trivial_filtering!(x::Optimizer, y::NodeBB)

    x._preprocess_termination_status = MOI.get(x.relaxed_optimizer, MOI.TerminationStatus())
    x._preprocess_result_status = MOI.get(x.relaxed_optimizer, MOI.PrimalStatus())
    valid_flag, feasible_flag = is_globally_optimal(x._preprocess_termination_status,
                                                    x._preprocess_result_status)

    if valid_flag
        if feasible_flag
            for j in 1:length(x._obbt_working_lower_index)
                @inbounds active_flag = x._obbt_working_lower_index[j]
                if active_flag
                    @inbounds vi = x._lower_variable_index[j]
                    diff = MOI.get(x.relaxed_optimizer, MOI.VariablePrimal(), vi)
                    @inbounds diff -= y.lower_variable_bounds[j]
                    if (abs(diff) <= x.obbt_tolerance)
                        @inbounds x._obbt_working_lower_index[j] = false
                    end
                end
            end
            for j in 1:length(x._obbt_working_upper_index)
                @inbounds active_flag = x._obbt_working_upper_index[j]
                if active_flag
                    @inbounds vi = x._lower_variable_index[j]
                    diff = MOI.get(x.relaxed_optimizer, MOI.VariablePrimal(), vi)
                    diff *= -1.0
                    @inbounds diff += y.upper_variable_bounds[j]
                    if (abs(diff) <= x.obbt_tolerance)
                        @inbounds x._obbt_working_upper_index[j] = false
                    end
                end
            end
        end
    end
    return
end

"""
$(FUNCTIONNAME)

Utility function used to set vector of booleans z to x & ~y. Avoids the
generation of conversion of the BitArray created by broadcasting logical operators.
"""
function bool_indx_diff(z::Vector{Bool},x::Vector{Bool}, y::Vector{Bool})
    for i = 1:length(z)
        @inbounds z[i] = (x[i] & ~y[i])
    end
    return
end

"""
$(FUNCTIONNAME)

Excludes OBBT on variable indices after a search in a filtering direction.
"""
function aggressive_filtering!(x::Optimizer, y::NodeBB)

    # Initial filtering vector (negative one direction per remark in Gleixner2017)
    variable_number = x._variable_number
    v = -ones(variable_number)

    # Copy prior index set (ignores linear and binary terms)
    obbt_var_len = length(x.obbt_variable_values)
    copyto!(x._old_low_index, x._obbt_working_lower_index)
    copyto!(x._old_upp_index, x._obbt_working_upper_index)
    copyto!(x._new_low_index, x._obbt_working_lower_index)
    copyto!(x._new_upp_index, x._obbt_working_upper_index)

    # Exclude unbounded directions
    for i in 1:obbt_var_len
        @inbounds active_flag = x._new_low_index[i]
        if active_flag
            @inbounds bnd = y.lower_variable_bounds[i]
            if (bnd == -Inf)
                @inbounds x._new_low_index[i] = false
            end
        end
    end
    for i in 1:obbt_var_len
        @inbounds active_flag = x._new_low_index[i]
        if active_flag
            @inbounds bnd = y.upper_variable_bounds[i]
            if bnd == Inf
                @inbounds x._new_low_index[i] = false
            end
        end
    end

    # Begin the main algorithm
    for k in 1:x.obbt_aggressive_max_iteration

        # Set index differences and vector for filtering direction
        bool_indx_diff(x._lower_indx_diff, x._old_low_index, x._new_low_index)
        bool_indx_diff(x._upper_indx_diff, x._old_upp_index, x._new_upp_index)

        for i in 1:obbt_var_len
            @inbounds active_flag = x._lower_indx_diff[i]
            if active_flag
                @inbounds bnd = v[i]
                if bnd < 0.0
                    @inbounds v[i] = 0.0
                end
            end
        end
        for i in 1:obbt_var_len
            @inbounds active_flag = x._upper_indx_diff[i]
            if active_flag
                @inbounds bnd = v[i]
                if (bnd > 0.0)
                    @inbounds v[i] = 0.0
                end
            end
        end

        # Termination Condition
        ((~any(x._new_low_index) & ~any(x._new_upp_index)) || (iszero(v))) && break
        if (k >= 2)
            if (count(x._lower_indx_diff) + count(x._upper_indx_diff)) < x.obbt_aggressive_min_dimension
                break
            end
        end

        # Set objective in OBBT problem to filtering vector
        MOI.set(x.relaxed_optimizer, MOI.ObjectiveSense(), MOI.MAX_SENSE)
        saf = SAF(SAT.(v, x._lower_variable_index), 0.0)
        MOI.set(x.relaxed_optimizer, MOI.ObjectiveFunction{SAF}(), saf)

        # Optimizes the problem and if successful filter additional bounds
        MOI.optimize!(x.relaxed_optimizer)

        x._preprocess_termination_status = MOI.get(x.relaxed_optimizer, MOI.TerminationStatus())
        x._preprocess_result_status = MOI.get(x.relaxed_optimizer, MOI.PrimalStatus())
        valid_flag, feasible_flag = is_globally_optimal(x._preprocess_termination_status,
                                                        x._preprocess_result_status)

        if valid_flag
            if feasible_flag
                variable_primal = MOI.get(x.relaxed_optimizer, MOI.VariablePrimal(), x._lower_variable_index)
                copyto!(x._new_low_index, x._old_low_index)
                copyto!(x._new_upp_index, x._old_upp_index)
                for i in 1:obbt_var_len
                    @inbounds active_flag = x._old_low_index[i]
                    if active_flag
                        @inbounds vp = variable_primal[i]
                        @inbounds xL = y.lower_variable_bounds[i]
                        if (vp == xL)
                            @inbounds x._new_low_index[i] = false
                        end
                    end
                end
                for i in 1:obbt_var_len
                    @inbounds active_flag = x._old_upp_index[i]
                    if active_flag
                        @inbounds vp = variable_primal[i]
                        @inbounds xU = y.upper_variable_bounds[i]
                        if (vp == xU)
                            @inbounds x._new_upp_index[i] = false
                        end
                    end
                end
            end
        else
            return false
        end
    end
    copyto!(x._obbt_working_lower_index, x._new_low_index)
    copyto!(x._obbt_working_upper_index, x._new_upp_index)
    return true
end

"""
$(FUNCTIONNAME)

Routine that determines if aggressive filtering should be used. Currently,
a user-specified option.
"""
aggressive_obbt_on_heurestic(x::Optimizer) = x.obbt_aggressive_on

"""
$(FUNCTIONNAME)

Performs OBBT with filtering and greedy ordering as detailed in:
Gleixner, A.M., Berthold, T., Müller, B. et al. J Glob Optim (2017) 67: 731.
https://doi.org/10.1007/s10898-016-0450-4
"""
function obbt(x::Optimizer)

    feasibility = true

    y = x._current_node
    x._current_xref .= 0.5*(y.upper_variable_bounds + y.lower_variable_bounds)

    # solve initial problem to feasibility
    update_relaxed_problem_box!(x, y)
    relax_problem!(x, x._current_xref, 1)
    relax_objective!(x, x._current_xref)
    objective_cut_linear!(x, 1)
    MOI.set(x.relaxed_optimizer, MOI.ObjectiveSense(), MOI.FEASIBILITY_SENSE)
    MOI.optimize!(x.relaxed_optimizer)

    # Sets indices to attempt OBBT on
    obbt_variables = x.obbt_variable_values
    copyto!(x._obbt_working_lower_index, obbt_variables)
    copyto!(x._obbt_working_upper_index, obbt_variables)

    # Prefiltering steps && and sets initial LP values
    trivial_filtering!(x, y)
    if aggressive_obbt_on_heurestic(x)
        feasibility = aggressive_filtering!(x, y)
    end

    x._preprocess_termination_status = MOI.get(x.relaxed_optimizer, MOI.TerminationStatus())
    x._preprocess_result_status = MOI.get(x.relaxed_optimizer, MOI.PrimalStatus())
    valid_flag, feasible_flag = is_globally_optimal(x._preprocess_termination_status,
                                                    x._preprocess_result_status)

    if valid_flag & feasible_flag
        xLP = MOI.get(x.relaxed_optimizer, MOI.VariablePrimal(), x._lower_variable_index)
    else
        return false
    end

    while (any(x._obbt_working_lower_index) || any(x._obbt_working_upper_index)) & ~isempty(y)

        # Get lower value
        lower_indx = -1
        upper_indx = -1
        lower_value = Inf
        upper_value = Inf

        # min of xLP - yL on active
        if any(x._obbt_working_lower_index)
            for i in 1:length(x._obbt_working_lower_index)
                @inbounds active_flag = x._obbt_working_lower_index[i]
                if active_flag
                    @inbounds temp_value = xLP[i]
                    @inbounds temp_value -= y.lower_variable_bounds[i]
                    if temp_value < lower_value
                        lower_value = temp_value
                        lower_indx = i
                    end
                end
            end
        end

        # min of yU - xLP on active
        if any(x._obbt_working_upper_index)
            for i in 1:length(x._obbt_working_upper_index)
                @inbounds active_flag = x._obbt_working_upper_index[i]
                if active_flag
                    @inbounds temp_value = y.upper_variable_bounds[i]
                    @inbounds temp_value -= xLP[i]
                    if temp_value < upper_value
                        upper_value = temp_value
                        upper_indx = i
                    end
                end
            end
        end

        if (lower_value <= upper_value)

            @inbounds x._obbt_working_lower_index[lower_indx] = false
            @inbounds var = x._lower_variable[lower_indx]
            MOI.set(x.relaxed_optimizer, MOI.ObjectiveSense(), MOI.MIN_SENSE)
            MOI.set(x.relaxed_optimizer, MOI.ObjectiveFunction{SV}(), var)
            MOI.optimize!(x.relaxed_optimizer)
            x._preprocess_termination_status = MOI.get(x.relaxed_optimizer, MOI.TerminationStatus())
            x._preprocess_result_status = MOI.get(x.relaxed_optimizer, MOI.PrimalStatus())
            valid_flag, feasible_flag = is_globally_optimal(x._preprocess_termination_status,
                                                            x._preprocess_result_status)

            if valid_flag
                if feasible_flag
                    xLP .= MOI.get(x.relaxed_optimizer, MOI.VariablePrimal(), x._lower_variable_index)
                    if is_integer_variable(x, lower_indx)
                        @inbounds y.lower_variable_bounds[lower_indx] = ceil(xLP[lower_indx])
                    else
                        @inbounds y.lower_variable_bounds[lower_indx] = xLP[lower_indx]
                    end
                    if isempty(y)
                        feasibility = false
                        break
                    end
                else
                    feasibility = false
                end
            else
                break
            end
        else

            x._obbt_working_upper_index[upper_indx] = false
            @inbounds var = x._lower_variable[upper_indx]
            MOI.set(x.relaxed_optimizer, MOI.ObjectiveSense(), MOI.MAX_SENSE)
            MOI.set(x.relaxed_optimizer, MOI.ObjectiveFunction{SV}(), var)
            MOI.optimize!(x.relaxed_optimizer)
            x._preprocess_termination_status = MOI.get(x.relaxed_optimizer, MOI.TerminationStatus())
            x._preprocess_result_status = MOI.get(x.relaxed_optimizer, MOI.PrimalStatus())
            valid_flag, feasible_flag = is_globally_optimal(x._preprocess_termination_status,
                                                            x._preprocess_result_status)

            if valid_flag
                if feasible_flag
                    xLP .= MOI.get(x.relaxed_optimizer, MOI.VariablePrimal(), x._lower_variable_index)
                    if is_integer_variable(x, upper_indx)
                        @inbounds y.upper_variable_bounds[upper_indx] = ceil(xLP[upper_indx])
                    else
                        @inbounds y.upper_variable_bounds[upper_indx] = xLP[upper_indx]
                    end
                    if isempty(y)
                        feasibility = false
                        break
                    end
                else
                    feasibility = false
                end
            else
                break
            end
        end
        trivial_filtering!(x, y)
    end

    return feasibility
end

"""
$(FUNCTIONNAME)

Performs the linear bound tightening.
"""
function lp_bound_tighten(m::Optimizer)

    feas = true
    n = m._current_node
    lvb = n.lower_variable_bounds
    uvb = n.upper_variable_bounds

    for i in 1:m.lp_repetitions

        # Runs Poor Man LP on constraints of form ax >= b
        for (func, constr, ind) in m._linear_geq_constraints
            if feas
                temp_value = -(constr.lower - func.constant)
                for term in func.terms
                    indx = term.variable_index.value
                    coeff = term.coefficient
                    @inbounds li = lvb[indx]
                    @inbounds ui = uvb[indx]
                    temp_value += max(coeff*ui, coeff*li)
                end
                for term in func.terms
                    if feas
                        indx = term.variable_index.value
                        coeff = term.coefficient
                        @inbounds li = lvb[indx]
                        @inbounds ui = uvb[indx]
                        term_value = max(coeff*ui, coeff*li)
                        cut_value = -(temp_value - term_value)/coeff
                        if (-coeff > 0.0)
                            if (li < cut_value)
                                if (cut_value < ui)
                                    @inbounds uvb[indx] = cut_value
                                end
                            else
                                feas = false
                                break
                            end
                        else
                            if (ui > cut_value)
                                if (cut_value > li)
                                    @inbounds lvb[indx] = cut_value
                                end
                            else
                                feas = false
                                break
                            end
                        end
                    end
                end
            else
                break
            end
        end

        # Runs Poor Man LP on constraints of form ax <= b
        for (func, constr, ind) in m._linear_leq_constraints
            if feas
                temp_value = (constr.upper - func.constant)
                for term in func.terms
                    indx = term.variable_index.value
                    coeff = term.coefficient
                    @inbounds li = lvb[indx]
                    @inbounds ui = uvb[indx]
                    temp_value -= min(coeff*ui, coeff*li)
                end
                for term in func.terms
                    if feas
                        indx = term.variable_index.value
                        coeff = term.coefficient
                        @inbounds li = lvb[indx]
                        @inbounds ui = uvb[indx]
                        term_value = min(coeff*ui, coeff*li)
                        cut_value = (temp_value + term_value)/coeff
                        if (coeff > 0.0)
                            if (li < cut_value)
                                if (cut_value < ui)
                                    @inbounds uvb[indx] = cut_value
                                end
                            else
                                feas = false
                                break
                            end
                        else
                            if (ui > cut_value)
                                if (cut_value > li)
                                    @inbounds lvb[indx] = cut_value
                                end
                            else
                                feas = false
                                break
                            end
                        end
                    end
                end
            else
                break
            end
        end

        for (func, constr, ind) in m._linear_eq_constraints
            if feas
                temp_value = (constr.value - func.constant)
                for term in func.terms
                    indx = term.variable_index.value
                    coeff = term.coefficient
                    @inbounds li = lvb[indx]
                    @inbounds ui = uvb[indx]
                    temp_value -= min(coeff*ui, coeff*li)
                end
                for term in func.terms
                    if feas
                        indx = term.variable_index.value
                        coeff = term.coefficient
                        @inbounds li = lvb[indx]
                        @inbounds ui = uvb[indx]
                        term_value = min(coeff*ui, coeff*li)
                        cut_value = (temp_value + term_value)/coeff
                        if (coeff > 0.0)
                            if (li < cut_value)
                                if (cut_value < ui)
                                    @inbounds uvb[indx] = cut_value
                                end
                            else
                                feas = false
                                break
                            end
                        else
                            if (ui > cut_value)
                                if (cut_value > li)
                                    @inbounds lvb[indx] = cut_value
                                end
                            else
                                feas = false
                                break
                            end
                        end
                    end
                end
            else
                break
            end

            if feas
                temp_value = -(constr.value - func.constant)
                for term in func.terms
                    indx = term.variable_index.value
                    coeff = term.coefficient
                    @inbounds li = lvb[indx]
                    @inbounds ui = uvb[indx]
                    temp_value += max(coeff*ui, coeff*li)
                end
                for term in func.terms
                    if feas
                        indx = term.variable_index.value
                        coeff = term.coefficient
                        @inbounds li = lvb[indx]
                        @inbounds ui = uvb[indx]
                        term_value = max(coeff*ui, coeff*li)
                        cut_value = -(temp_value - term_value)/coeff
                        if (-coeff > 0.0)
                            if (li < cut_value)
                                if (cut_value < ui)
                                    @inbounds uvb[indx] = cut_value
                                end
                            else
                                feas = false
                                break
                            end
                        else
                            if (ui > cut_value)
                                if (cut_value > li)
                                    @inbounds lvb[indx] = cut_value
                                end
                            else
                                feas = false
                                break
                            end
                        end
                    end
                end
            else
                break
            end
        end

        (~feas) && (break)
    end

    m._current_node = NodeBB(lvb, uvb, n.lower_bound, n.upper_bound, n.depth, n.id)

    return feas
end

function check_univariate_quad(f::MOI.ScalarQuadraticFunction{Float64})
    (length(f.affine_terms) == 1) &&
    (length(f.quadratic_terms) == 1) &&
    (f.affine_terms[1].variable_index == f.quadratic_terms[1].variable_index_1 == f.quadratic_terms[1].variable_index_2)
end

get_value(set::MOI.LessThan{Float64}) = set.upper
get_value(set::MOI.GreaterThan{Float64}) = set.lower
get_value(set::MOI.EqualTo{Float64}) = set.value

function get_univariate_coeff(func::MOI.ScalarQuadraticFunction{Float64}, set::T) where {T<:MOI.AbstractScalarSet}
    a = func.quadratic_terms[1].coefficient
    b = (length(func.affine_terms) > 0) ?  func.affine_terms[1].coefficient : 0.0
    c = get_value(set) - func.constant
    vi = func.quadratic_terms[1].variable_index_1.value
    a,b,c,vi
end

# Checks to see if constraint is a bivariant quadratic term
function check_bivariate_quad(f::MOI.ScalarQuadraticFunction{Float64})
    vIndx = Int64[]
    (length(f.quadratic_terms) > 3) && (return false)
    (length(f.affine_terms) > 2) && (return false)
    for i in f.affine_terms push!(vIndx,i.variable_index.value) end
    for i in f.quadratic_terms push!(vIndx,i.variable_index_1.value, i.variable_index_2.value) end
    vIndx = Int64[]
    unique_vIndx = unique(vIndx)
    unique_cnt = length(unique_vIndx)
    if (unique_cnt == 2)
        return true, unique_vIndx[1], unique_vIndx[2]
    else
        return false, nothing, nothing
    end
end

function get_bivariate_coeff(func::MOI.ScalarQuadraticFunction{Float64},set::T,vxvalue::Int,vyvalue::Int) where {T<:MOI.AbstractScalarSet}
    acnt = length(func.affine_terms)
    (vxvalue != nothing) && (vx = MOI.VariableIndex(vxvalue))
    (vyvalue != nothing) && (vy = MOI.VariableIndex(vyvalue))

    for qd_term in func.quadratic_terms
        if (qd_term.variable_index1 == vx && qd_term.variable_index2 == vx)
            ax = qd_term.coefficient
        elseif (qd_term.variable_index1 == vy && qd_term.variable_index2 == vy)
            ay = qd_term.coefficient
        else
            axy = qd_term.coefficient
        end
    end

    affine_coefficient_1 = func.affine_terms[1].coefficient
    affine_coefficient_2 = func.affine_terms[2].coefficient
    if acnt == 2
        if (func.affine_terms[1].variable_index1 == vx)
            bx = affine_coefficient_1
            by = affine_coefficient_2
        else
            bx = affine_coefficient_2
            by = affine_coefficient_1
        end
    else
        if (func.affine_terms[1].variable_index1 == vx)
            bx = affine_coefficient_1
            by = 0.0
        else
            bx = 0.0
            by = affine_coefficient_1
        end
    end
    c = get_value(set) - func.constant
end

"""
$(FUNCTIONNAME)

Classifies constraints as univariate or bivariate and adds
them to storage vectors.
"""
function classify_quadratics!(m::Optimizer)
    a = 0.0
    b = 0.0
    c = 0.0
    # Check for Univariate and Bivariate Lesser Constraints
    for (func,set,indx) in m._quadratic_leq_constraints
        if check_univariate_quad(func)
            a,b,c,vi = get_univariate_coeff(func,set)
            a_neg = -1.0*a
            b_neg = -1.0*b
            c_neg = -1.0*c
            push!(m._univariate_quadratic_leq_constraints,(a_neg,b_neg,c_neg,vi))
        else
             flag,vxi,vyi = check_bivariate_quad(func)
             if flag
                ax,ay,axy,bx,by,c = get_bivariate_coeff(func,set)
                push!(m._bivariate_quadratic_geq_constraints,(-ax,-ay,-axy,-bx,-by,c,vxi,vyi))
            end
        end
    end

    # Check for Univariate and Bivariate Greater Constraints
    for (func,set,indx) in m._quadratic_geq_constraints
        if check_univariate_quad(func)
            a,b,c,vi = get_univariate_coeff(func,set)
            push!(m._univariate_quadratic_geq_constraints,(a,b,c,vi))
        else
            flag,vxi,vyi = check_bivariate_quad(func)
            if flag
               ax,ay,axy,bx,by,c = get_bivariate_coeff(func,set)
               push!(m._bivariate_quadratic_geq_constraints,(ax,ay,axy,bx,by,-c,vxi,vyi))
           end
        end
    end

    # Check for Univariate and Bivariate Equality Constraints
    for (func,set,indx) in m._quadratic_eq_constraints
        if check_univariate_quad(func)
            a,b,c,vi = get_univariate_coeff(func,set)
            push!(m._univariate_quadratic_eq_constraints,(a,b,c,vi))
        else
            flag,vxi,vyi = check_bivariate_quad(func)
            if flag
               ax,ay,axy,bx,by,c = get_bivariate_coeff(func,set)
               push!(m._bivariate_quadratic_geq_constraints,(ax,ay,axy,bx,by,-c,vxi,vyi))
           end
           flag,vxi,vyi = check_bivariate_quad(func)
           if flag
              ax,ay,axy,bx,by,c = get_bivariate_coeff(func,set)
              push!(m._bivariate_quadratic_geq_constraints,(-ax,-ay,-axy,-bx,-by,c,vxi,vyi))
          end
        end
    end
    return
end

"""
$(FUNCTIONNAME)

Kernel of the bound tightening operation on univariant qudaratic functions.
Called for each univariate function.
"""
function univariate_kernel(m::Optimizer,a::Float64,b::Float64,c::Float64,vi::Int)
        flag = true
        term1 = c + (b^2)/(4.0*a)
        term2 = term1/a
        if ((term1 > 0.0) && (a < 0.0)) # No solution, fathom node
            flag = false
        elseif (term2 >= 0.0)
            xlo = m._current_node.lower_variable_bounds[vi]
            xhi = m._current_node.upper_variable_bounds[vi]
            chk1 = -sqrt(term2)-b/(2.0*a)
            chk2 = sqrt(term2)-b/(2.0*a)
            if (a > 0.0)
                (chk1 < xlo) && (m._current_node.lower_variable_bounds[vi] = max(xlo,chk2))
                (chk2 > xhi) && (m._current_node.upper_variable_bounds[vi] = min(xhi,chk1))
            else
                m._current_node.lower_variable_bounds[vi] = max(xlo,chk1)
                m._current_node.upper_variable_bounds[vi] = min(xhi,chk2)
            end
            if (m._current_node.lower_variable_bounds[vi] <= m._current_node.upper_variable_bounds[vi])
                flag = true
            end
        else
            flag = true
        end
        return flag
end

"""
$(FUNCTIONNAME)

Performs bound tightening on all univariate quadratic functions.
"""
function univariate_quadratic(m::Optimizer)
    feas = true
    # fathom ax^2 + bx + c > l quadratics
    for (a, b, c, vi) in m._univariate_quadratic_geq_constraints
        if ~feas
            break
        end
        feas = univariate_kernel(m, a, b, c , vi)
    end
    # fathom ax^2 + bx + c < u quadratics
    for (a, b, c, vi) in m._univariate_quadratic_leq_constraints
        if ~feas
            break
        end
        feas = univariate_kernel(m, -a, -b, -c, vi)
    end
    # fathom ax^2 + bx + c = v quadratics
    for (a, b, c, vi) in m._univariate_quadratic_eq_constraints
        if ~feas
            break
        end
        feas = univariate_kernel(m, a, b, c, vi)
        if ~feas
            break
        end
        feas = univariate_kernel(m, -a, -b, -c, vi)
     end
     return feas
end

#=
"""
$(FUNCTIONNAME)

Kernel of the bound tightening operation on bivariate qudaratic functions.
Called for each bivariate function.
"""
function bivariate_kernel(m::Optimizer,n::NodeBB,ax::Float64,ay::Float64,axy::Float64,
                         bx::Float64,by::Float64,vi1::Int,vi2)
        # Case distinction from Vigerske disseration (TO DO)
end

"""
$(FUNCTIONNAME)

Performs bound tightening on all bivariate quadratic functions.
"""
function bivariate_quadratic(m::Optimizer,n::NodeBB)
    feas = true
    # fathom ax^2 + bx + c > l quadratics
    for (ax,ay,axy,bx,by,vi) in m.bivariate_quadratic_geq_constraints
        feas = bivariate_kernel(m,n,ax,ay,axy,bx,by,vi)
        (~feas) && return feas
    end
    return feas
end
=#

"""
$(FUNCTIONNAME)

Performs forward-reverse pass on directed graph as
part of constraint propagation.
"""
function cpwalk(x::Optimizer)

    n = x._current_node
    evaluator = x._relaxed_evaluator

    # runs at midpoint bound
    midx = (n.upper_variable_bounds + n.lower_variable_bounds)/2.0

    # set working node to n, copies pass parameters from EAGO optimizer
    evaluator.current_node = n
    evaluator.has_reverse = true
    prior_sg_tighten = evaluator.subgrad_tighten

    evaluator.subgrad_tighten = false
    evaluator.cp_repetitions = x.cp_repetitions
    evaluator.cp_tolerance = x.cp_tolerance

    # Run forward-reverse pass and retreive node for interval forward-reverse pass
    evaluator.subgrad_tighten = ~x.cp_interval_only
    feas = forward_reverse_pass(evaluator, midx)
    @inbounds n.lower_variable_bounds[:] = evaluator.current_node.lower_variable_bounds
    @inbounds n.upper_variable_bounds[:] = evaluator.current_node.upper_variable_bounds

    # resets forward reverse scheme for lower bounding problem
    evaluator.has_reverse = false
    evaluator.subgrad_tighten = prior_sg_tighten

    return feas
end
