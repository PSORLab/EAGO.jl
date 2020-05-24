function buffer_kernel(func::SQF)
    buffer = OrderedDict{Int, Float64}()
    for term in func.quadratic_terms
        buffer[term.variable_index_1.value] = 0.0
        buffer[term.variable_index_2.value] = 0.0
    end
    for term in func.affine_terms
        buffer[term.variable_index.value] = 0.0
    end
    saf1 = SAF([SAT(0.0, VI(k)) for k in keys(buffer)], 0.0)
    saf2 = SAF([SAT(0.0, VI(k)) for k in keys(buffer)], 0.0)
    nx = length(buffer)
    buffer, saf1, saf2, nx
end

function add_buffered_quadratic!(m::Optimizer, func::SQF, set::LT)
    buffer, saf1, saf2, nx = buffer_kernel(func)
    cfunc = copy(func)
    cfunc.constant -= set.upper
    push!(m._inner_constraints._sqf_leq, BufferQuadratic(cfunc, QD_LT, buffer, saf1, saf2, nx))
end

function add_buffered_quadratic!(m::Optimizer, func::SQF, set::GT)
    buffer, saf1, saf2, nx = buffer_kernel(func)
    cfunc = copy(func)
    cfunc.constant = set.lower - cfunc.constant
    for i = 1:length(cfunc.quadratic_terms)
        term = @inbounds cfunc.quadratic_terms[i]
        @inbounds cfunc.quadratic_terms[i] = SQT(-term.coefficient, term.variable_index_1, term.variable_index_2)
    end
    for i = 1:length(cfunc.affine_terms)
        term = @inbounds cfunc.affine_terms[i]
        @inbounds cfunc.affine_terms[i] = SAT(-term.coefficient, term.variable_index)
    end
    push!(m._inner_constraints._sqf_leq, BufferQuadratic(cfunc, QD_GT, buffer, saf1, saf2, nx))
end

function add_buffered_quadratic!(m::Optimizer, func::SQF, set::ET, indx::Int)
    buffer_a, saf1_a, saf2_a, nx_a = buffer_kernel(func)
    buffer_b, saf1_b, saf2_b, nx_b = buffer_kernel(func)
    cfunc1 = copy(func)
    cfunc1.constant -= set.upper
    cfunc2 = copy(func)
    cfunc2.constant = set.value - cfunc.constant
    for i = 1:length(cfunc2.quadratic_terms)
        term = @inbounds cfunc2.quadratic_terms[i]
        @inbounds cfunc2.quadratic_terms[i] = SQT(-term.coefficient, term.variable_index_1, term.variable_index_2)
    end
    for i = 1:length(cfunc2.affine_terms)
        term = @inbounds cfunc2.affine_terms[i]
        @inbounds cfunc2.affine_terms[i] = SAT(-term.coefficient, term.variable_index)
    end
    push!(m._inner_constraints._sqf_leq, BufferQuadratic(cfunc1, QD_LT, buffer_a, saf1_a, saf2_a, nx_a))
    push!(m._inner_constraints._sqf_leq, BufferQuadratic(cfunc2, QD_GT, buffer_b, saf1_b, saf2_b, nx_b))
end
