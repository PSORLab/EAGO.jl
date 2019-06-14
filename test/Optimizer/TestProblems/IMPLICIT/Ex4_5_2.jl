#=
@testset "Begin Implicit Solver 4.5.2 Example" begin
end
=#

opt = EAGO.Optimizer(absolute_tolerance = 0.001, relative_tolerance = 0.001)

function obj(z,p)

    alpha = 2.5; beta = 0.52;
    c_steam = 21.67*10^(-3); c_cool = 4.65*10^(-3);

    y_4A = one(p[1]) - z[8] - z[9]

    C_1_cap = 132718 + z[3]*(369*z[4] - 1113.9*y[5])
    C_2_cap = 25000 + z[7]*(6984.5*z[8]-3869.53*z[9]^2)

    C_CSTR_cap = 25764 + 8178*p[1]

    C_1_op = z[3]*(3.0 + 3.611*y_4A + 7.71*z[8])*(c_steam + c_cool)
    C_2_op = z[4]*(26.21 + 29.45*z[8])*(c_steam + c_cool)

    C_ann = (C_1_cap + C_2_cap + C_CSTR_cap)/alpha + beta*(C_1_op + C_2_op)

    return C_ann
end

# z = (F_1, - > 1
#      F_2,  - > 2
#      F_3,  - > 3
#      y_3A,  - > 4
#      y_3B,  - > 5
#      y_3C,  - > 6
#      F_4,  - > 7
#      y_4B,  - > 8
#      y_4c, - > 9
#      F_6,  - > 10
#      F_7 - > 11)
function h!(out,z,p)

    v_A = 8.937*10^(-2); v_B = 1.018*10^(-1); v_C = 1.130*10^(-1);
    F5 = 50;

    r_1 = p[2]*z[4]/(z[4]*v_A + z[5]*v_B + z[6]*v_C);
    r_2 = p[3]*z[5]/(z[4]*v_A + z[5]*v_B + z[6]*v_C);

    out[1] = z[1] + z[11] - z[2]
    out[2] = z[2] - z[4]*z[3] - p[1]*r_1
    out[3] = p[1]*(r_1 - r_2) - F5
    out[4] = p[1]*r_2 - z[6]*z[3]
    out[5] = one(p[1]) - z[4] - z[5] - z[6]
    out[6] = z[3] - z[7] - z[11]
    out[7] = z[4]*z[3] - z[11]
    out[8] = z[5]*z[3] - z[8]*z[7]
    out[9] = z[7] - F5 - z[10]
    out[10] = z[8]*z[7] - F5
    out[11] = z[9]*z[7] - z[10]
end

# TO DO:
function hjac!(out,x,p)

    v_A = 8.937*10^(-2); v_B = 1.018*10^(-1); v_C = 1.130*10^(-1);

    dr_1_yA =
    dr_1_yB =
    dr_1_yC =

    dr_2_yA =
    dr_2_yB =
    dr_2_yC =

    out[1,1] = one(p[1])
    out[1,2] = -one(p[1])
    out[1,11] = one(p[1])

    out[2,2] =
    out[2,3] =
    out[2,4] =

    out[3] =
    out[4] =
    out[5] =
    out[6] =
    out[7] =
    out[8] =
    out[9] =
    out[10] =
    out[11] =

end

pl = [14.0 0.405 0.0545]; pu = [24.0 0.415 0.0555]
xl = [50.0 100.0 100.0 0.5 0.1 0.001 50.0 0.142857 0.142857*10^(-3) 0.0 40.0]
xu = [150.0 300.0 300.0 0.9 0.5 0.2 60.0 1.0 1.0 10.0 250.0]
var, opt = solve_implicit(obj, h!, xl, xu, pl, pu, opt = opt, hj = hjac!)

pval = JuMP.value.(var)
fval = JuMP.objective_value(opt)
tstatus = JuMP.termination_status(opt)
pstatus = JuMP.primal_status(opt)

p_bool1 = isapprox(pval, 21.68, atol = 0.01)
p_bool2 = isapprox(pval, 0.415, atol = 0.01)
p_bool3 = isapprox(pval, 0.0545, atol = 0.01)
f_bool = isapprox(fval, 289780, atol = 0.001)
t_bool = (tstatus == MOI.OPTIMAL)
p_bool = (pstatus == MOI.FEASIBLE_POINT)
