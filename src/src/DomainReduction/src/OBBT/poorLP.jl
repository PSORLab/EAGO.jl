# Linear Program Contractor
function poorLP_contractor(f,g,X::Array{Interval},UBD,nx::Int64)
    # gets prior McCormick relaxation setting
    mu_temp = EAGOSmoothMcCormickGrad.MC_param.mu
    set_diff_relax!(0)

    # sets variable bounds
    l = [X[i].lo for i=1:nx]
    u = [X[i].lo for i=1:nx]
    l_out = copy(l)
    u_out = copy(u)

    # computes relaxation of g and f
    x_mc = [SMCg{nx,Float64}(mid.(X[i]),mid.(X[i]),seed(nx,i),seed(nx,i),X[i],false,X,mid.(X))]
    f_mc = f(x)
    g_mc = g(x)

    ### constructs linear problem
    A = zeros(length(g_mc)+1,nx)
    b = zeros(length(g_mc)+1)
    A[1,:] = [f.cv_grad[i] for i=1:nx]
    b[1] = UBD - f.cv + sum([f.cv_grad[i]*mid.(X[i]) for i=1:nx])
    A[2:length(g_mc)+1,:] = [g_mc[i].cv_grad[j] for i=1:length(g_mc),j=1:nx]
    b[2:length(g_mc)+1] = [sum([g_mc[i].cv_grad[j]*mid.(X[i]) for j=1:nx])-g_mc[i].cv for i=1:length(g_mc)]

    for hind=1:nx
        for i=1:(length(g_mc)+1)
            if (A[i,hind] > 0.0)
                upper_bnd = (1.0/A[i,hind])*(b[i]-sum([ hind == j ? 0.0 : min(A[i,j]*u_out[j],A[i,j]*l_out[j]) for j=1:nx]))
                if (u_out[hind] > upper_bnd)
                    u_out[hind] = upper_bnd
                end
            elseif (A[i,hind] < 0.0)
                lower_bnd = (1.0/A[i,hind])*(b[i]-sum([ hind == j ? 0.0 : min(A[i,j]*u_out[j],A[i,j]*l_out[j]) for j=1:nx]))
                if (l_out[hind] < lower_bnd)
                    u_out[hind] = upper_bnd
                end
            end
        end
        if (lout[hind]>uout[hind])
            return X, false
        end
    end

    # reset McCormick relaxation to original value
    set_diff_relax!(mu_temp)

    temp_arr = [Interval(l_out[i],u_out[i]) for i=1:nx]
    return temp_arr,true
end
