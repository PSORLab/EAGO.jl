"""
    is_convex_quadratic

Returns true if `func` < 0  based on eigenvalue tests, false otherwise.
"""
function is_convex_quadratic(func::MOI.ScalarQuadraticFunction{Float64},NumVar::Int,mult::Float64)
    # Riguous Convexity Test
    Q = spzeros(NumVar,NumVar)
    for term in func.quadratic_terms
        if term.coefficient != 0.0
            Q[term.variable_index_1.value,term.variable_index_2.value] = mult*term.coefficient
        end
    end
    if length(Q.nzval) > 1
        eigmin = LinearAlgebra.eigmin(Array(Q))
        if (eigmin) > 0.0
            return true
        end
    else
        if Q.nzval[1] > 0.0
            return true
        else
            return false
        end
    end
    return false
end


"""
    quadratic_convexity!

Assigns boolean value to constraint_convexity dictionary entry corresponding to
constraint index that is true if constraint is shown to be convex and false
otherwise.
"""
function quadratic_convexity!(src::Optimizer)
    for (func, set, ind) in src.quadratic_leq_constraints
        MOIindx = MOI.ConstraintIndex{MOI.ScalarQuadraticFunction{Float64},MOI.LessThan{Float64}}(ind)
        src.constraint_convexity[MOIindx] = is_convex_quadratic(func,src.variable_number,1.0)
    end
    for (func,set,ind) in src.quadratic_geq_constraints
        MOIindx = MOI.ConstraintIndex{MOI.ScalarQuadraticFunction{Float64},MOI.GreaterThan{Float64}}(ind)
        src.constraint_convexity[MOIindx] = is_convex_quadratic(func,src.variable_number,-1.0)
    end
    for (func,set,ind) in src.quadratic_eq_constraints
        MOIindx = MOI.ConstraintIndex{MOI.ScalarQuadraticFunction{Float64},MOI.EqualTo{Float64}}(ind)
        src.constraint_convexity[MOIindx] = false
    end
end

"""
    relax_convex_kernel
"""
function relax_convex_kernel(func::MOI.ScalarQuadraticFunction{Float64}, n::NodeBB, src::Optimizer, x0::Vector{Float64})
    VarNum = length(n)
    temp_sto = zeros(Float64,VarNum)
    terms_coeff = Float64[]
    terms_index = Int[]
    term_constant = func.constant
    for term in func.quadratic_terms
        temp_constant = x0[src.variable_index_to_storage[term.variable_index_1.value]]
        temp_constant *= x0[src.variable_index_to_storage[term.variable_index_2.value]]
        term_constant -= term.coefficient*temp_constant
    end

    # Adds affine terms
    for term in func.affine_terms
        temp_sto[src.variable_index_to_storage[term.variable_index.value]] += term.coefficient
    end

    # Adds quadratic terms
    for term in func.quadratic_terms
        @inbounds temp_sto[src.variable_index_to_storage[term.variable_index_1.value]] += 2.0*term.coefficient*x0[src.variable_index_to_storage[term.variable_index_1.value]]
    end

    temp = 0.0
    for (i,term) in enumerate(temp_sto)
        if (term != 0.0)
            push!(terms_coeff,term)
            push!(terms_index,src.storage_index_to_variable[i])
        end
    end

    varIndx = [MOI.VariableIndex(src.variable_index_to_storage[i]) for i in terms_index]

    return varIndx, terms_coeff, term_constant
end
"""
    relax_convex_quadratic_inner!

Default routine for relaxing convex quadratic constraint `lower` < `func` < `upper`
on node `n`. Takes affine bounds of convex part at point `x0`.
"""
function relax_convex_quadratic_inner!(trg, src::Optimizer, func::MOI.ScalarQuadraticFunction{Float64},
                                 lower::Float64, upper::Float64, n::NodeBB, x0::Vector{Float64})


    varIndx, terms_coeff, term_constant = relax_convex_kernel(func, n, src, x0)

    if (lower == upper)
        saf_term = MOI.ScalarAffineTerm{Float64}.(terms_coeff,varIndx)
        saf = MOI.ScalarAffineFunction{Float64}(saf_term,term_constant)
        MOI.add_constraint(trg,saf,MOI.LessThan{Float64}(upper))
        term_constant *= -1.0
        for i in 1:length(terms_coeff)
            @inbounds terms_coeff[i] *= -1.0
        end
        saf1_term = MOI.ScalarAffineTerm{Float64}.(terms_coeff,varIndx)
        saf1 = MOI.ScalarAffineFunction{Float64}(saf1_term,term_constant)
        MOI.add_constraint(trg,saf1,MOI.LessThan{Float64}(-lower))
    elseif (lower != -Inf)
        term_constant *= -1.0
        for i in 1:length(terms_coeff)
            @inbounds terms_coeff[i] *= -1.0
        end
        saf_term = MOI.ScalarAffineTerm{Float64}.(terms_coeff,varIndx)
        saf = MOI.ScalarAffineFunction{Float64}(saf_term,term_constant)
        MOI.add_constraint(trg,saf,MOI.LessThan{Float64}(-lower))
    else
        saf_term = MOI.ScalarAffineTerm{Float64}.(terms_coeff,varIndx)
        saf = MOI.ScalarAffineFunction{Float64}(saf_term,term_constant)
        MOI.add_constraint(trg,saf,MOI.LessThan{Float64}(upper))
    end
end

"""
    relax_nonconvex_kernel

Default calculation kernel for nonconvex relaxation.
"""
function relax_nonconvex_kernel(func,n,src,x0)
    terms_coeff = Float64[]
    terms_index = Int[]
    VarNum = length(n)
    quadratic_coefficient = zeros(Float64,VarNum)
    quadratic_constant = func.constant
    for term in func.quadratic_terms
        # TODO: Move a, NodeIndx1, NodeIndx2 to separate storage
        a = term.coefficient
        NodeIndx1 = src.variable_index_to_storage[term.variable_index_1.value]
        NodeIndx2 = src.variable_index_to_storage[term.variable_index_2.value]
        xL1 = n.lower_variable_bounds[NodeIndx1]
        xU1 = n.upper_variable_bounds[NodeIndx1]
        x01 = x0[NodeIndx1]
        if NodeIndx1 == NodeIndx2               # quadratic terms
            if (a > 0.0)
                quadratic_coefficient[NodeIndx1] = 2.0*a*x01
                quadratic_constant -= 2.0*x01*x01
            else
                quadratic_coefficient[NodeIndx1] = a*(xL1 + xU1)
                quadratic_constant -= a*xL1*xU1
            end
        else
            xL2 = n.lower_variable_bounds[NodeIndx2]
            xU2 = n.upper_variable_bounds[NodeIndx2]
            if (a > 0.0)
                check = (xU1 - xL1)*x0[NodeIndx2] + (xU2 - xL2)*x0[NodeIndx1] - xU1*xU2 + xL1*xL2 <= 0.0
                if check
                    quadratic_coefficient[NodeIndx1] += a*xL2
                    quadratic_coefficient[NodeIndx2] += a*xL1
                    quadratic_constant -= a*xL1*xL2
                else
                    quadratic_coefficient[NodeIndx1] += a*xU2
                    quadratic_coefficient[NodeIndx2] += a*xU1
                    quadratic_constant -= a*xU1*xU2
                end
            else
                check = (xU1 - xL1)*x0[NodeIndx2] + (xU2 - xL2)*x0[NodeIndx1] <= xU1*xL2 - xL1*xU2
                if check
                    quadratic_coefficient[NodeIndx1] += a*xL2
                    quadratic_coefficient[NodeIndx2] += a*xU1
                    quadratic_constant -= a*xU1*xL2
                else
                    quadratic_coefficient[NodeIndx1] += a*xU2
                    quadratic_coefficient[NodeIndx2] += a*xL1
                    quadratic_constant -= a*xL1*xU2
                end
            end
        end
    end
    # Adds affine terms
    # TODO: Move term in func.affine_terms to separate storage
    for term in func.affine_terms
        NodeIndx3 = src.variable_index_to_storage[term.variable_index.value]
        quadratic_coefficient[NodeIndx3] += term.coefficient
    end

    temp = 0.0
    for (i,term) in enumerate(quadratic_coefficient)
        if (term != 0.0)
            push!(terms_coeff,term)
            push!(terms_index,src.storage_index_to_variable[i])
        end
    end

    varIndx = [MOI.VariableIndex(src.variable_index_to_storage[i]) for i in terms_index]

    return varIndx, terms_coeff, quadratic_constant
end

"""
    relax_nonconvex_quadratic!

Default routine for relaxing nonconvex quadratic constraint `lower` < `func` < `upper`
on node `n`. Takes affine bounds of convex part at point `x0` and secant line
bounds on concave parts.
"""
function relax_nonconvex_quadratic!(trg, src::Optimizer, func::MOI.ScalarQuadraticFunction{Float64},
                                    lower::Float64, upper::Float64, n::NodeBB, x0::Vector{Float64})

    varIndx, terms_coeff, quadratic_constant = relax_nonconvex_kernel(func,n,src,x0)

    if (lower == upper)
        saf_term = MOI.ScalarAffineTerm{Float64}.(terms_coeff,varIndx)
        saf = MOI.ScalarAffineFunction{Float64}(saf_term,0.0)
        MOI.add_constraint(trg,saf,MOI.LessThan{Float64}(upper-quadratic_constant))
        varIndx, terms_coeff, quadratic_constant = relax_nonconvex_kernel(func,n,src,x0)
        quadratic_constant *= -1.0
        for i in 1:length(terms_coeff)
            terms_coeff[i] *= -1.0
        end
        saf1_term = MOI.ScalarAffineTerm{Float64}.(terms_coeff,varIndx)
        saf1 = MOI.ScalarAffineFunction{Float64}(saf1_term,0.0)
        set1 = MOI.LessThan{Float64}(-lower-quadratic_constant)
        MOI.add_constraint(trg,saf1,set1)
    elseif (lower != -Inf)
        quadratic_constant *= -1.0
        for i in 1:length(terms_coeff)
            terms_coeff[i] *= -1.0
        end
        saf_term = MOI.ScalarAffineTerm{Float64}.(terms_coeff,varIndx)
        saf = MOI.ScalarAffineFunction{Float64}(saf_term,0.0)
        MOI.add_constraint(trg,saf,MOI.LessThan{Float64}(-lower-quadratic_constant))
    else
        saf_term = MOI.ScalarAffineTerm{Float64}.(terms_coeff,varIndx)
        saf = MOI.ScalarAffineFunction{Float64}(saf_term,0.0)
        MOI.add_constraint(trg,saf,MOI.LessThan{Float64}(upper-quadratic_constant))
    end
end

"""
    relax_quadratic!

Relaxes all quadratic constraints in `src` optimizer over the domain of `n` using
the relaxation scheme `r` if it is defined or by forming affine bounds and adding
them to the `trg` optimizer.
"""
function relax_quadratic!(trg, src::Optimizer, n::NodeBB, r::RelaxationScheme)

    x0 = (n.upper_variable_bounds - n.lower_variable_bounds)/2.0

    # Relax Convex Constraint Terms
    for (func, set, ind) in src.quadratic_leq_constraints
        MOIindx = MOI.ConstraintIndex{MOI.ScalarQuadraticFunction{Float64},MOI.LessThan{Float64}}(ind)
        if src.constraint_convexity[MOIindx]
            relax_convex_quadratic_inner!(trg, src, func, -Inf, set.upper, n, x0)
        else
            relax_nonconvex_quadratic!(trg, src, func, -Inf, set.upper, n, x0)
        end
    end
    for (func,set,ind) in src.quadratic_geq_constraints
        MOIindx = MOI.ConstraintIndex{MOI.ScalarQuadraticFunction{Float64},MOI.GreaterThan{Float64}}(ind)
        if src.constraint_convexity[MOIindx]
            relax_convex_quadratic_inner!(trg, src, func, set.lower, Inf, n, x0)
        else
            relax_nonconvex_quadratic!(trg, src, func, set.lower, Inf, n, x0)
        end
    end
    for (func,set,ind) in src.quadratic_eq_constraints
        MOIindx = MOI.ConstraintIndex{MOI.ScalarQuadraticFunction{Float64},MOI.EqualTo{Float64}}(ind)
        if src.constraint_convexity[MOIindx]
            relax_convex_quadratic_inner!(trg, src, func, set.value, set.value, n, x0)
        else
            relax_nonconvex_quadratic!(trg, src, func, set.value, set.value, n, x0)
        end
    end

    # Relax quadratic objective
    if isa(src.objective, MOI.ScalarQuadraticFunction{Float64}) # quadratic objective
        if (src.optimization_sense == MOI.MAX_SENSE)
            objf = src.objective
            m_objf= MOI.ScalarQuadraticFunction{Float64}(objf.affine_terms, objf.quadratic_terms, -objf.constant)
            for i in 1:length(m_objf.affine_terms)
                m_objf.affine_terms[i] = MOI.ScalarAffineTerm{Float64}(-m_objf.affine_terms[i].coefficient, m_objf.affine_terms[i].variable_index)
            end
            for i in 1:length(m_objf.quadratic_terms)
                m_objf.quadratic_terms[i] =  MOI.ScalarQuadraticTerm{Float64}(-m_objf.quadratic_terms[i].coefficient, m_objf.quadratic_terms[i].variable_index_1, m_objf.quadratic_terms[i].variable_index_2)
            end
            varIndx, terms_coeff, quadratic_constant = relax_nonconvex_kernel(m_objf, n, src, x0)
            saf_term = MOI.ScalarAffineTerm{Float64}.(terms_coeff, varIndx)
            obj_func = MOI.ScalarAffineFunction{Float64}(saf_term, quadratic_constant)
            MOI.set(trg, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(), obj_func)
        else
            if (src.objective_convexity)
                varIndx, terms_coeff, quadratic_constant = relax_convex_kernel(func, n, src, x0)
            else
                varIndx, terms_coeff, quadratic_constant = relax_nonconvex_kernel(src.objective, n, src, x0)
            end
            saf_term = MOI.ScalarAffineTerm{Float64}.(terms_coeff, varIndx)
            obj_func = MOI.ScalarAffineFunction{Float64}(saf_term, quadratic_constant)
            MOI.set(trg, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(), obj_func)
        end
    end
    #MOI.set(trg, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),  QuadMidPoint)
end
