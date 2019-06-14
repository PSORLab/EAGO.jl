function univariate_kernel(n::NodeBB,a::Float64,b::Float64,c::Float64,vi::Int)
        term1 = c + (b^2)/(4.0*a)
        term2 = term1/a
        if ((term1 > 0.0) && (a < 0.0)) # No solution, fathom node
            return false
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
                return true
            end
        else
            return true
        end
end

function univariate_quadratic(m::Optimizer,n::NodeBB)
    feas = true
    # fathom ax^2 + bx + c > l quadratics
    for (a,b,c,vi) in m.univariate_quadratic_geq_constraints
        feas = univariate_kernel(n,a,b,c,vi)
        (feas == false) && return feas
    end
    # fathom ax^2 + bx + c < u quadratics
    for (a,b,c,vi) in m.univariate_quadratic_leq_constraints
         feas = univariate_kernel(n,-a,-b,-c,vi)
         (feas == false) && return feas
    end
    # fathom ax^2 + bx + c = v quadratics
    for (a,b,c,vi) in m.univariate_quadratic_eq_constraints
          feas = univariate_kernel(n,a,b,c,vi)
          (feas == false) && return feas
          feas = univariate_kernel(n,-a,-b,-c,vi)
          (feas == false) && return feas
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

```
Checks to see if constraint is a univariant quadratic term
```
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

    # Should be done
    for qd_term in func.quadratic_terms
        if (qd_term.variable_index1 == vx && qd_term.variable_index2 == vx)
            ax = qd_term.coefficient
        elseif (qd_term.variable_index1 == vy && qd_term.variable_index2 == vy)
            ay = qd_term.coefficient
        else
            axy = qd_term.coefficient
        end
    end

    # assigns affine terms... (DONE)
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

```
Classifies constraints as univariate or bivariate and adds them to storage vector
```
function classify_quadratics!(m::Optimizer)
    a = 0.0
    b = 0.0
    c = 0.0
    # Check for Univariate and Bivariate Lesser Constraints
    for (func,set,indx) in m.quadratic_leq_constraints
        if check_univariate_quad(func)
            a,b,c,vi = get_univariate_coeff(func,set)
            a_neg = -1.0*a
            b_neg = -1.0*b
            c_neg = -1.0*c
            push!(m.univariate_quadratic_leq_constraints,(a_neg,b_neg,c_neg,vi))
        else
             flag,vxi,vyi = check_bivariate_quad(func)
             if flag
                ax,ay,axy,bx,by,c = get_bivariate_coeff(func,set)
                push!(m.bivariate_quadratic_geq_constraints,(-ax,-ay,-axy,-bx,-by,c,vxi,vyi))
            end
        end
    end

    # Check for Univariate and Bivariate Greater Constraints
    for (func,set,indx) in m.quadratic_geq_constraints
        if check_univariate_quad(func)
            a,b,c,vi = get_univariate_coeff(func,set)
            push!(m.univariate_quadratic_geq_constraints,(a,b,c,vi))
        else
            flag,vxi,vyi = check_bivariate_quad(func)
            if flag
               ax,ay,axy,bx,by,c = get_bivariate_coeff(func,set)
               push!(m.bivariate_quadratic_geq_constraints,(ax,ay,axy,bx,by,-c,vxi,vyi))
           end
        end
    end

    # Check for Univariate and Bivariate Equality Constraints
    for (func,set,indx) in m.quadratic_eq_constraints
        if check_univariate_quad(func)
            a,b,c,vi = get_univariate_coeff(func,set)
            push!(m.univariate_quadratic_eq_constraints,(a,b,c,vi))
        else
            flag,vxi,vyi = check_bivariate_quad(func)
            if flag
               ax,ay,axy,bx,by,c = get_bivariate_coeff(func,set)
               push!(m.bivariate_quadratic_geq_constraints,(ax,ay,axy,bx,by,-c,vxi,vyi))
           end
           flag,vxi,vyi = check_bivariate_quad(func)
           if flag
              ax,ay,axy,bx,by,c = get_bivariate_coeff(func,set)
              push!(m.bivariate_quadratic_geq_constraints,(-ax,-ay,-axy,-bx,-by,c,vxi,vyi))
          end
        end
    end
end
