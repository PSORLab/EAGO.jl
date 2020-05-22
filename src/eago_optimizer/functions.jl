mutable struct BufferedQuadratic{S <: INEQ_SETS}
    func::SQF
    set::S
    buffer::OrderedDict{Int,Float64}
    saf1::SAF
    saf2::SAF
    nx::Int
    ci::CI{SQF, S}
end
function BufferedQuadratic(func::SQF, set::S, ci::CI{SQF,S}) where {S <: INEQ_SETS}
    buffer = OrderedDict{Int, Float64}()
    for term in func.quadratic_terms
        buffer[term.variable_index_1.value] = 0.0
        buffer[term.variable_index_2.value] = 0.0
    end
    for term in func.affine_terms
        buffer[term.variable_index.value] = 0.0
    end
    for term in buffer
        saf1 = SAF([SAT(0.0, VI(k)) for k in keys(buffer)], 0.0)
        saf2 = SAF([SAT(0.0, VI(k)) for k in keys(buffer)], 0.0)
    end
    nx = length(buffer)
    ci = CI{SQF,S}(0)
    BufferQuadratic(func, set, buffer, saf1, saf2, nx, ci)
end

function affine_relax_quadratic!(func::SQF, buffer::OrderedDict{Int,Float64}, saf::SAF, n::NodeBB, x::Vector{Float64}, p1::Bool)

    lower_bounds = n.lower_variable_bounds
    upper_bounds = n.upper_variable_bounds

    quadratic_constant = func.constant

    for term in func.quadratic_terms

        a = p1 ? term.coefficient : -term.coefficient
        idx1 = term.variable_index_1.value
        idx2 = term.variable_index_2.value

        x0_1 = @inbounds x0[idx1]
        xL_1 = @inbounds lower_bounds[idx1]
        xU_1 = @inbounds upper_bounds[idx1]

        if idx1 === idx2

            @inbounds buffer[idx1] = (a > 0.0) ? 2.0*a*x0_1 : a*(xL_1 + xU_1)
            quadratic_constant -= (a > 0.0) ? x0_1*x0_1 : a*xL_1*xU_1

        else
            x0_2 = @inbounds x0[idx2]
            xL_2 = @inbounds lower_bounds[idx2]
            xU_2 = @inbounds upper_bounds[idx2]

            if a > 0.0
                check_ref = (xU_1 - xL_1)*x0_2 + (xU_2 - xL_2)*x0_1
                if check_ref <= xU_1*xU_2 - xL_1*xL_2
                    @inbounds buffer[idx1] += a*xL_2
                    @inbounds buffer[idx2] += a*xL_1
                    quadratic_constant -= a*xL_1*xL_2
                else
                    @inbounds buffer[idx1] += a*xU_2
                    @inbounds buffer[idx2] += a*xU_1
                    quadratic_constant -= a*xU_1*xU_2
                end
            else
                check_ref = (xU_1 - xL_1)*x0_2 - (xU_2 - xL_2)*x0_1
                if check_ref <= xU_1*xL_2 - xL_1*xU_2
                    @inbounds buffer[idx1] += a*xL_2
                    @inbounds buffer[idx2] += a*xU_1
                    quadratic_constant -= a*xU_1*xL_2
                else
                    @inbounds buffer[idx1] += a*xU_2
                    @inbounds buffer[idx2] += a*xL_1
                    quadratic_constant -= a*xL_1*xU_2
                end
            end
        end
    end

    for term in func.affine_terms
        a = p1 ? term.coefficient : -term.coefficient
        idx = term.variable_index.value
        @inbounds buffer[idx] += a
    end

    count = 1
    for (key, value) in buffer
        @inbounds saf.affine_terms[count] = SAT(value, VI(key))
        count += 1
    end
    saf.constant = quadratic_constant

    return
end

affine_relax!(f::BufferedQuadratic{LT}, n::NodeBB, x::Vector{Float64}) = affine_relax_quadratic!(f.func, f.buffer, f.saf1, n, x, true)
affine_relax!(f::BufferedQuadratic{GT}, n::NodeBB, x::Vector{Float64}) = affine_relax_quadratic!(f.func, f.buffer, f.saf1, n, x, false)
function affine_relax!(f::BufferedQuadratic{ET}, n::NodeBB, x::Vector{Float64})
    affine_relax_quadratic!(f.func, f.buffer, f.saf1, n, x, true)
    affine_relax_quadratic!(f.func, f.buffer, f.saf2, n, x, false)
    nothing
end

relax!(m::Optimizer, f::BufferedQuadratic, k::Int)
