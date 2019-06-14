#=
@testset "Begin Implicit Solver 5.2 Example" begin
end
=#

opt = EAGO.Optimizer(absolute_tolerance = 0.001, relative_tolerance = 0.001)
a = [37.3692 18.5805 6.25]
c = [0.602 1.211 3.6]

function obj(z,p)
    out = zero(p[1])
    temp = zero(p[1])
    for j in 1:3
        temp = a[j]*(p[j]-c[j])^2
        for i in 1:3
            if (i != j)
                temp += a[i]*(p[i]-c[i])
            end
        end
        temp_inner = zero(p[1])
        for i in 1:3
            temp_inner += z[i]*(-1)^(i+1)
        end
        temp -= 5.0*((j-1)*(j-2)*(z[2]-z[1]) + temp_inner)
        out += temp^2
    end
    return out
end

function h!(out,z,p)

    out[1] = (1.00*10^(-9))*(exp(38.0*z[1])-1.0) + p[1]*z[1] - 1.6722*z[2] + 0.6689*z[3] - 8.0267
    out[2] = (1.98*10^(-9))*(exp(38.0*z[2])-1.0) + 0.6622*z[1] + p[2]*z[2] + 0.6622*z[3] + 4.0535
    out[3]  = (1.00*10^(-9))*(exp(38.0*z[3])-1.0) + z[1] - z[2] + p[3]*z[3] - 6.0
end

function hjac!(out,x,p)

     out[1,1] = 38.0*(1.00*10^(-9))*exp(38.0*z[1]) + p[1]
     out[1,2] = 1.6722*one(p[1])
     out[1,3] = 0.6689*one(p[1])

     out[2,1] = 0.6622*one(p[1])
     out[2,2] = 38.0*(1.98*10^(-9))*exp(38.0*z[2]) + p[2]
     out[2,3] = 0.6622*one(p[1])

     out[3,1] = one(p[1])
     out[3,2] = -one(p[1])
     out[3,3] = 38.0*(1.00*10^(-9))*exp(38.0*z[3])+ p[3]

end

pl = [0.6020 1.2110 3.6]; pu = [0.7358 1.4801 4.4];
xl = [0.5180 -3.9748 0.3296]; xu = [0.5847 -3.0464 0.5827]
var, opt = solve_implicit(obj, h!, xl, xu, pl, pu, opt = opt, hj = hjac!)

pval = JuMP.value.(var)
fval = JuMP.objective_value(opt)
tstatus = JuMP.termination_status(opt)
pstatus = JuMP.primal_status(opt)

p_bool1 = isapprox(pval, 0.703918, atol = 0.01)
p_bool2 = isapprox(pval, 1.43648, atol = 0.01)
p_bool3 = isapprox(pval, 3.61133, atol = 0.01)
f_bool = isapprox(fval, 626.565, atol = 0.001)
t_bool = (tstatus == MOI.OPTIMAL)
p_bool = (pstatus == MOI.FEASIBLE_POINT)
