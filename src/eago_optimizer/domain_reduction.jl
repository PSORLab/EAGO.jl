"""
    variable_dbbt!

Tightens the bounds of the `_current_node` using the current global upper bound
and the duality information obtained from the relaxation.
"""
function variable_dbbt!(x::NodeBB, mult_lo::Vector{Float64}, mult_hi::Vector{Float64},
                        LBD::Float64, UBD::Float64, nx::Int)

    cut = 0.0
    vb = 0.0
    if LBD <= UBD
        for i = 1:nx
            @inbounds mult_lo = mult_lo[i]
            @inbounds mult_hi = mult_hi[i]
            if mult_lo > 0.0
                @inbounds cut = x.upper_variable_bounds[i] - (UBD - LBD)/mult_lo
                @inbounds vb = x.lower_variable_bounds[i]
                if cut > vb
                    @inbounds x.lower_variable_bounds[i] = cut
                end
             elseif mult_hi > 0.0
                 @inbounds cut = x.lower_variable_bounds[i] + (UBD - LBD)/mult_hi
                 @inbounds vb = x.upper_variable_bounds[i]
                 if cut < vb
                     @inbounds x.upper_variable_bounds[i] = cut
                 end
             end
         end
    end

    return
end

"""
    trivial_filtering!

Excludes OBBT on variable indices that are tight for the solution of the relaxation.
"""
function trival_filtering!(x::Optimizer, y::NodeBB,
                           obbt_working_lower_index::Vector{Int64},
                           obbt_working_upper_index::Vector{Int64})

    termination_status = MOI.get(x.relaxed_optimizer, MOI.TerminationStatus())
    result_status_code = MOI.get(x.relaxed_optimizer, MOI.PrimalStatus())
    valid_flag, feasible_flag = is_globally_optimal(termination_status, result_status_code)

    if valid_flag
        if feasible_flag
            low_delete = Int64[]
            upp_delete = Int64[]
            for i in obbt_working_lower_index
                @inbounds vi = x._lower_variable_index[i]
                diff = MOI.get(x.relaxed_optimizer, MOI.VariablePrimal(), vi)
                @inbounds diff -= y.lower_variable_bounds[i]
                (abs(diff) <= x.obbt_tolerance) && push!(low_delete, i)
            end
            for i in obbt_working_upper_index
                @inbounds vi = x._lower_variable_index[i]
                diff = MOI.get(x.relaxed_optimizer, MOI.VariablePrimal(), vi)
                diff *= -1.0
                @inbounds diff += y.upper_variable_bounds[i]
                (abs(diff) <= x.obbt_tolerance) && push!(upp_delete, i)
            end
            deleteat!(obbt_working_lower_index, low_delete)
            deleteat!(obbt_working_upper_index, upp_delete)
        end
    end
    return
end

"""
    aggressive_filtering!

Excludes OBBT on variable indices after a search in a filtering direction.
"""
function aggressive_filtering!(x::Optimizer, y::NodeBB,
                               obbt_working_lower_index::Vector{Int64},
                               obbt_working_upper_index::Vector{Int64})

    # Initial filtering vector (negative one direction per remark in Gleixner2017)
    variable_number = x._variable_number
    v = -ones(variable_number)

    # Copy prior index set (ignores linear and binary terms)
    old_low_index = copy(obbt_working_lower_index)
    old_upp_index = copy(obbt_working_upper_index)
    new_low_index = copy(obbt_working_lower_index)
    new_upp_index = copy(obbt_working_upper_index)

    # Exclude unbounded directions
    low_delete = Int64[]
    upp_delete = Int64[]
    for i in new_low_index
        @inbounds bnd = y.lower_variable_bounds[i]
        if ~(bnd == -Inf)
            push!(low_delete, i)
        end
    end
    for i in new_upp_index
        @inbounds bnd = y.upper_variable_bounds[i]
        if ~(bnd == Inf)
            push!(upp_delete, i)
        end
    end
    deleteat!(new_low_index, low_delete)
    deleteat!(new_upp_index, upp_delete)
    empty!(low_delete)
    empty!(upp_delete)

    # Begin the main algorithm
    for k in 1:x.obbt_aggressive_max_iteration

        # Set index differences and vector for filtering direction
        lower_indx_diff = setdiff(old_low_index, new_low_index)
        upper_indx_diff = setdiff(old_upp_index, new_upp_index)
        low_delete = Int64[]
        upp_delete = Int64[]

        for i in lower_indx_diff
            @inbounds bnd = v[i]
            if bnd < 0.0
                @inbounds v[i] = 0.0
            end
        end
        for i in upper_indx_diff
            @inbounds bnd = v[i]
            if (bnd > 0.0)
                @inbounds v[i] = 0.0
            end
        end

        # Termination Condition
        ((isempty(new_low_index) & isempty(new_upp_index)) || (iszero(v))) && break
        if (k >= 2)
            if (length(lower_indx_diff) + length(upper_indx_diff)) < x.obbt_aggressive_min_dimension
                break
            end
        end

        # Set objective in OBBT problem to filtering vector
        MOI.set(x.relaxed_optimizer, MOI.ObjectiveSense(), MOI.MAX_SENSE)
        saf = SAF(SAT.(v, x._lower_variable_index), 0.0)
        MOI.set(x.relaxed_optimizer, MOI.ObjectiveFunction{SAF}(), saf)

        # Optimizes the problem and if successful filter additional bounds
        MOI.optimize!(x.relaxed_optimizer)

        termination_status = MOI.get(x.relaxed_optimizer, MOI.TerminationStatus())
        result_status_code = MOI.get(x.relaxed_optimizer, MOI.PrimalStatus())
        valid_flag, feasible_flag = is_globally_optimal(termination_status, result_status_code)

        if valid_flag
            if feasible_flag
                variable_primal = MOI.get(x.relaxed_optimizer, MOI.VariablePrimal(), x._lower_variable_index)
                new_low_index = copy(old_low_index)
                new_upp_index = copy(old_upp_index)
                for i in old_low_index
                    @inbounds chk = (variable_primal[i] == y.lower_variable_bounds[i])
                    chk && push!(low_delete, i)
                end
                for i in old_upp_index
                    @inbounds (variable_primal[i] == y.upper_variable_bounds[i])
                    chk && push!(upp_delete, i)
                end
                deleteat!(new_low_index, low_delete)
                deleteat!(new_upp_index, upp_delete)
            end
        else
            return false
        end
    end
    obbt_working_lower_index = new_low_index
    obbt_working_upper_index = new_upp_index
    return true
end
aggressive_obbt_on_heurestic(x::Optimizer) = x.obbt_aggressive_on

"""
    obbt

Performs OBBT with filtering and greedy ordering
"""
function obbt(x::Optimizer)

    feasibility = true

    y = x._current_node
    ymid = 0.5*(y.upper_variable_bounds + y.lower_variable_bounds)

    # solve initial problem to feasibility
    update_relaxed_problem_box!(x, y)
    relax_problem!(x.ext_type, x, ymid)
    relax_objective!(x.ext_type, x, ymid)
    objective_cut_linear!(x)
    MOI.set(x.relaxed_optimizer, MOI.ObjectiveSense(), MOI.FEASIBILITY_SENSE)
    MOI.optimize!(x.relaxed_optimizer)

    # Sets indices to attempt OBBT on
    obbt_variables = x._lower_variable_values
    obbt_working_lower_index = Int64[]
    obbt_working_upper_index = Int64[]
    append!(obbt_working_lower_index, obbt_variables)
    append!(obbt_working_upper_index, obbt_variables)

    # Prefiltering steps && and sets initial LP values
    trival_filtering!(x, y, obbt_working_lower_index,
                            obbt_working_upper_index)

    if aggressive_obbt_on_heurestic(x)
        feasibility = aggressive_filtering!(x, y, obbt_working_lower_index,
                                                  obbt_working_upper_index)
    end
    xLP = MOI.get(x.relaxed_optimizer, MOI.VariablePrimal(), x._lower_variable_index)

    while ~(isempty(obbt_working_lower_index) && isempty(obbt_working_upper_index)) && ~isempty(y)

        # Get lower value
        lower_indx = 0
        upper_indx = 0
        if (isempty(obbt_working_lower_index))
            lower_value = Inf
        else
            @inbounds find_indx = obbt_working_lower_index[1]
            @inbounds lower_value = xLP[find_indx]
            @inbounds lower_value -= y.lower_variable_bounds[find_indx]
            lower_indx = 1
            for i in 2:length(obbt_working_lower_index)
                @inbounds find_indx = obbt_working_lower_index[i]
                @inbounds temp_value = xLP[find_indx]
                @inbounds temp_value -= y.lower_variable_bounds[find_indx]
                if temp_value < lower_value
                    lower_value = temp_value
                    lower_indx = i
                end
            end
        end
        if (isempty(obbt_working_upper_index))
            upper_value = Inf
        else
            @inbounds find_indx = obbt_working_upper_index[1]
            @inbounds upper_value = y.upper_variable_bounds[find_indx]
            @inbounds upper_value -= xLP[find_indx]
            upper_indx = 1
            for i in 2:length(obbt_working_upper_index)
                @inbounds find_indx = obbt_working_upper_index[i]
                @inbounds temp_value = y.upper_variable_bounds[find_indx]
                @inbounds temp_value -= xLP[find_indx]
                if temp_value < upper_value
                    upper_value = temp_value
                    upper_indx = i
                end
            end
        end

        if (lower_value <= upper_value)

            deleteat!(obbt_working_lower_index, lower_indx)

            # set objectives
            @inbounds var = x._lower_variable[lower_indx]
            MOI.set(x.relaxed_optimizer, MOI.ObjectiveSense(), MOI.MIN_SENSE)
            MOI.set(x.relaxed_optimizer, MOI.ObjectiveFunction{SV}(), var)
            MOI.optimize!(x.relaxed_optimizer)
            termination_status = MOI.get(x.relaxed_optimizer, MOI.TerminationStatus())
            result_status_code = MOI.get(x.relaxed_optimizer, MOI.PrimalStatus())
            valid_flag, feasible_flag = is_globally_optimal(termination_status, result_status_code)

            if valid_flag
                if feasible_flag
                    @inbounds xLP[:] = MOI.get(x.relaxed_optimizer, MOI.VariablePrimal(), x._lower_variable_index)
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

            # Get index
            deleteat!(obbt_working_upper_index, upper_indx)

            @inbounds var = x._lower_variable[upper_indx]
            MOI.set(x.relaxed_optimizer, MOI.ObjectiveSense(), MOI.MAX_SENSE)
            MOI.set(x.relaxed_optimizer, MOI.ObjectiveFunction{SV}(), var)
            MOI.optimize!(x.relaxed_optimizer)
            termination_status = MOI.get(x.relaxed_optimizer, MOI.TerminationStatus())
            result_status_code = MOI.get(x.relaxed_optimizer, MOI.PrimalStatus())
            valid_flag, feasible_flag = is_globally_optimal(termination_status, result_status_code)

            if valid_flag
                if feasible_flag
                    @inbounds xLP[:] = MOI.get(x.relaxed_optimizer, MOI.VariablePrimal(), x._lower_variable_index)
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
        #GenerateLVB!(x)
        trival_filtering!(x, y, obbt_working_lower_index,
                                obbt_working_upper_index)
    end

    return feasibility
end

"""
    lp_bound_tighten

Performs the linear bound tightening.
"""
function lp_bound_tighten(m::Optimizer)

    n = m._current_node
    lvb = n.lower_variable_bounds
    uvb = n.upper_variable_bounds

    # Runs Poor Man LP on constraints of form ax >= b
    for (func, constr, ind) in m._linear_geq_constraints
        temp_value = (constr.lower - func.constant)
        for term in func.terms
            indx = term.variable_index.value
            coeff = term.coefficient
            @inbounds li = lvb[indx]
            @inbounds ui = uvb[indx]
            temp_value += -max(coeff*ui, coeff*li)
        end
        for term in func.terms
            indx = term.variable_index.value
            coeff = term.coefficient
            @inbounds li = lvb[indx]
            @inbounds ui = uvb[indx]
            term_value = -max(coeff*ui, coeff*li)
            cut_value = (temp_value - term_value)/coeff
            if (coeff > 0.0)
                if (li < cut_value)
                    (cut_value > ui) && (return false)
                    @inbounds n.lower_variable_bounds[vi] = cut_value
                end
            else
                if (ui > cut_value)
                    (cut_value < li) && (return false)
                    @inbounds n.upper_variable_bounds[indx] = cut_value
                end
            end
        end
    end

    # Runs Poor Man LP on constraints of form ax <= b
    for (func, constr, ind) in m._linear_leq_constraints
        temp_value = (constr.upper - func.constant)
        for term in func.terms
            indx = term.variable_index.value
            coeff = term.coefficient
            @inbounds li = lvb[indx]
            @inbounds ui = uvb[indx]
            temp_value += -min(coeff*ui, coeff*li)
        end
        for term in func.terms
            indx = term.variable_index.value
            coeff = term.coefficient
            @inbounds li = lvb[indx]
            @inbounds ui = uvb[indx]
            term_value = -min(coeff*ui, coeff*li)
            cut_value = (temp_value - term_value)/coeff
            if (term.coefficient < 0.0 )
                if (li < cut_value)
                    (cut_value > ui) && (return false)
                    @inbounds n.lower_variable_bounds[vi] = cut_value
                end
            else
                if (ui > cut_value)
                    (cut_value < li) && (return false)
                    @inbounds n.upper_variable_bounds[vi] = cut_value
                end
            end
        end
    end

    for (func, constr, ind) in m._linear_eq_constraints
        temp_value = (constr.value - func.constant)
        for term in func.terms
            indx = term.variable_index.value
            coeff = term.coefficient
            @inbounds li = lvb[indx]
            @inbounds ui = uvb[indx]
            temp_value += -max(coeff*ui, coeff*li)
        end
        for term in func.terms
            indx = term.variable_index.value
            coeff = term.coefficient
            @inbounds li = lvb[indx]
            @inbounds ui = uvb[indx]
            term_value = -max(coeff*ui, coeff*li)
            cut_value = (temp_value - term_value)/coeff
            if (coeff > 0.0 )
                if (li < cut_value)
                    (cut_value > ui) && (return false)
                    @inbounds n.lower_variable_bounds[vi] = cut_value
                end
            else
                if (ui > cut_value)
                    (cut_value < li) && (return false)
                    @inbounds n.upper_variable_bounds[vi] = cut_value
                end
            end
        end
        temp_value = (constr.value - func.constant)
        for term in func.terms
            indx = term.variable_index.value
            coeff = term.coefficient
            @inbounds li = lvb[indx]
            @inbounds ui = uvb[indx]
            TempValue += -min(coeff*ui, coeff*li)
        end
        for term in func.terms
            indx = term.variable_index.value
            coeff = term.coefficient
            @inbounds li = lvb[indx]
            @inbounds ui = uvb[indx]
            term_value = -min(coeff*ui, coeff*li)
            cut_value = (temp_value - term_value)/coeff
            if (coeff < 0.0 )
                if (li < cut_value)
                    (cut_value > ui) && (return false)
                    @inbounds n.lower_variable_bounds[vi] = cut_value
                end
            else
                if (ui > cut_value)
                    (cut_value < li) && (return false)
                    @inbounds n.upper_variable_bounds[vi] = cut_value
                end
            end
        end
    end

    return true
end

function check_univariate_quad(f::MOI.ScalarQuadraticFunction{Float64})
    (length(f.affine_terms) == 1) &&
    (length(f.quadratic_terms) == 1) &&
    (f.affine_terms[1].variable_index == f.quadratic_terms[1].variable_index_1 == f.quadratic_terms[1].variable_index_2)
end

get_value(set::MOI.LessThan{Float64}) = set.upper
get_value(set::MOI.GreaterThan{Float64}) = set.lower
get_value(set::MOI.EqualTo{Float64}) = set.value

function get_univariate_coeff(func::MOI.ScalarQuadraticFunction{Float64},set::T) where {T<:MOI.AbstractScalarSet}
    a = func.quadratic_terms[1].coefficient
    b = (length(func.affine_terms) > 0) ?  func.affine_terms[1].coefficient : 0.0
    c = get_value(set) - func.constant
    vi = func.quadratic_terms[1].variable_index_1.value
    a,b,c,vi
end

```
Checks to see if constraint is a bivariant quadratic term
```
function check_bivariate_quad(f::MOI.ScalarQuadraticFunction{Float64})
    vIndx = Int[]
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
    classify_quadratics!

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

function univariate_kernel(n::NodeBB,a::Float64,b::Float64,c::Float64,vi::Int)
        flag = true
        term1 = c + (b^2)/(4.0*a)
        term2 = term1/a
        if ((term1 > 0.0) && (a < 0.0)) # No solution, fathom node
            flag = false
        elseif (term2 >= 0.0)
            xlo = n.lower_variable_bounds[vi]
            xhi = n.upper_variable_bounds[vi]
            chk1 = -sqrt(term2)-b/(2.0*a)
            chk2 = sqrt(term2)-b/(2.0*a)
            if (a > 0.0)
                (chk1 < xlo) && (n.lower_variable_bounds[vi] = max(xlo,chk2))
                (chk2 > xhi) && (n.upper_variable_bounds[vi] = min(xhi,chk1))
            else
                n.lower_variable_bounds[vi] = max(xlo,chk1)
                n.upper_variable_bounds[vi] = min(xhi,chk2)
            end
            if (n.lower_variable_bounds[vi] <= n.upper_variable_bounds[vi])
                flag = true
            end
        else
            flag = true
        end
        return flag
end

function univariate_quadratic(m::Optimizer)
    feas = true
    # fathom ax^2 + bx + c > l quadratics
    for (a, b, c, vi) in m._univariate_quadratic_geq_constraints
        if ~feas
            break
        end
        feas = univariate_kernel(n, a, b, c , vi)
    end
    # fathom ax^2 + bx + c < u quadratics
    for (a, b, c, vi) in m._univariate_quadratic_leq_constraints
        if ~feas
            break
        end
        feas = univariate_kernel(n, -a, -b, -c, vi)
    end
    # fathom ax^2 + bx + c = v quadratics
    for (a, b, c, vi) in m._univariate_quadratic_eq_constraints
        if ~feas
            break
        end
        feas = univariate_kernel(n, a, b, c, vi)
        if ~feas
            break
        end
        feas = univariate_kernel(n, -a, -b, -c, vi)
     end
     return feas
end

function bivariate_kernel(m::Optimizer,n::NodeBB,ax::Float64,ay::Float64,axy::Float64,
                         bx::Float64,by::Float64,vi1::Int,vi2)
        # Case distinction from Vigerske disseration (TO DO)
end

function bivariate_quadratic(m::Optimizer,n::NodeBB)
    feas = true
    # fathom ax^2 + bx + c > l quadratics
    for (ax,ay,axy,bx,by,vi) in m.bivariate_quadratic_geq_constraints
        feas = bivariate_kernel(m,n,ax,ay,axy,bx,by,vi)
        (~feas) && return feas
    end
    return feas
end

"""
    cpwalk

Performs forward-reverse pass on directed graph as
part of constraint propagation.
"""
function cpwalk(x::Optimizer)

    n = x._current_node
    evaluator = x.working_evaluator_block.evaluator

    # runs at midpoint bound
    midx = (n.upper_variable_bounds + n.lower_variable_bounds)/2.0

    # set working node to n, copies pass parameters from EAGO optimizer
    evaluator.current_node = n
    evaluator.has_reverse = true
    prior_sg_tighten = evaluator.subgrad_tighten

    evaluator.subgrad_tighten = false
    evaluator.cp_reptitions = x.cp_interval_reptitions
    evaluator.cp_tolerance = x.cp_interval_tolerance

    # Run forward-reverse pass and retreive node for interval forward-reverse pass
    feas = forward_reverse_pass(evaluator, midx)
    @inbounds n.lower_variable_bounds[:] = evaluator.current_node.lower_variable_bounds
    @inbounds n.upper_variable_bounds[:] = evaluator.current_node.upper_variable_bounds

    # Run forward-reverse pass and retreive node for mccormick forward-reverse pass
    if feas
        evaluator.subgrad_tighten = true
        evaluator.cp_reptitions = x.cp_mccormick_reptitions
        evaluator.cp_tolerance = x.cp_mccormick_tolerance

        feas = forward_reverse_pass(evaluator, midx)
        @inbounds n.lower_variable_bounds[:] = evaluator.current_node.lower_variable_bounds
        @inbounds n.upper_variable_bounds[:] = evaluator.current_node.upper_variable_bounds
    end

    # resets forward reverse scheme for lower bounding problem
    evaluator.has_reverse = false
    evaluator.subgrad_tighten = prior_sg_tighten

    return feas
end
