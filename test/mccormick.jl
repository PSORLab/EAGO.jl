# Create functions for comparing MC object to reference object and detail failures
function chk_ref_kernel(y::MC{N,T}, yref::MC{N,T}, mctol::Float64) where {N, T<:EAGO.RelaxTag}
    pass_flag::Bool = true
    descr = "Failing Components: ("
    ~isapprox(y.cv, yref.cv; atol = mctol) && (descr = descr*" cv = $(y.cv)"; pass_flag = false)
    ~isapprox(y.cc, yref.cc; atol = mctol) && (descr = descr*" cc = $(y.cc) "; pass_flag = false)
    ~isapprox(y.Intv.lo, yref.Intv.lo; atol = mctol) && (descr = descr*" Intv.lo = $(y.Intv.lo) "; pass_flag = false)
    ~isapprox(y.Intv.hi, yref.Intv.hi; atol = mctol) && (descr = descr*" Intv.hi = $(y.Intv.hi) "; pass_flag = false)
    ~isapprox(y.cv_grad[1], yref.cv_grad[1]; atol = mctol) && (descr = descr*" cv_grad[1] = $(y.cv_grad[1]) "; pass_flag = false)
    ~isapprox(y.cv_grad[2], yref.cv_grad[2]; atol = mctol) && (descr = descr*" cv_grad[2] = $(y.cv_grad[2]) "; pass_flag = false)
    ~isapprox(y.cc_grad[1], yref.cc_grad[1]; atol = mctol) && (descr = descr*" cc_grad[1] = $(y.cc_grad[1]) "; pass_flag = false)
    ~isapprox(y.cc_grad[2], yref.cc_grad[2]; atol = mctol) && (descr = descr*" cc_grad[2] = $(y.cc_grad[2]) "; pass_flag = false)
    (descr !== "Failing Components: (") && println(descr*")")
    pass_flag
end
check_vs_ref1(f::Function, x::MC, yref::MC, mctol::Float64) = chk_ref_kernel(f(x), yref, mctol)
check_vs_ref2(f::Function, x::MC, y::MC, yref::MC, mctol::Float64) = chk_ref_kernel(f(x,y), yref, mctol)
check_vs_refv(f::Function, x::MC, c::Float64, yref::MC, mctol::Float64) = chk_ref_kernel(f(x, c), yref, mctol)
check_vs_refv(f::Function, x::MC, c::Int, yref::MC, mctol::Float64) = chk_ref_kernel(f(x, c), yref, mctol)
function check_vs_ref2(f::Function, x::MC, c::Float64, yref1::MC, yref2::MC, mctol::Float64)
    pass_flag = check_vs_ref_kernel(f(x, c), yref1, mctol)
    check_vs_ref_kernel(f(c, x), yref2, mctol) && pass_flag
end


@testset "Test Univariate" begin

   mctol = 1E-4

   ##### Exponent #####
   x_exp = MC{2,NS}(2.0, 3.0, Interval{Float64}(1.0,4.0), seed_gradient(1,Val(2)), seed_gradient(1,Val(2)), false)
   x_exp2 = MC{2,NS}(2.0, 3.0, Interval{Float64}(1.0,4.0), seed_gradient(1,Val(2)), seed_gradient(1,Val(2)), false)
   x_exp10 = MC{2,NS}(2.0, 3.0, Interval{Float64}(1.0,4.0), seed_gradient(1,Val(2)), seed_gradient(1,Val(2)), false)
   x_expm1 = MC{2,NS}(2.0, 3.0, Interval{Float64}(1.0,4.0), seed_gradient(1,Val(2)), seed_gradient(1,Val(2)), false)

   yref_exp = MC{2,NS}(7.38906, 37.3049, Interval{Float64}(2.71828, 54.5982), SVector{2,Float64}([7.38906, 0.0]), SVector{2,Float64}([17.2933, 0.0]), false)
   yref_exp2 = MC{2,NS}(4.0, 11.33333333333333, Interval{Float64}(1.99999, 16.0001), SVector{2,Float64}([2.77259, 0.0]), SVector{2,Float64}([4.66667, 0.0]), false)
   yref_exp10 = MC{2,NS}(100.0, 6670.0000000000055, Interval{Float64}(9.999999999999999999, 10000.00000000001), SVector{2,Float64}([230.25850929940458, 0.0]), SVector{2,Float64}([3330.0, 0.0]), false)
   yref_expm1 = MC{2,NS}(6.38905609893065, 36.304860631582514, Interval{Float64}(1.71828, 53.5982), SVector{2,Float64}([7.38906, 0.0]), SVector{2,Float64}([17.2933, 0.0]), false)

   @test check_vs_ref1(exp, x_exp, yref_exp, mctol)
   @test check_vs_ref1(exp2, x_exp2, yref_exp2, mctol)
   @test check_vs_ref1(exp10, x_exp10, yref_exp10, mctol)
   @test check_vs_ref1(expm1, x_expm1, yref_expm1, mctol)

   ##### Logarithm #####
   x_log = MC{2,NS}(2.0, 3.0, Interval{Float64}(1.0,4.0), seed_gradient(1,Val(2)), seed_gradient(1,Val(2)), false)
   x_log2 = MC{2,NS}(2.0, 3.0, Interval{Float64}(1.0,4.0), seed_gradient(1,Val(2)), seed_gradient(1,Val(2)), false)
   x_log10 = MC{2,NS}(2.0, 3.0, Interval{Float64}(1.0,4.0), seed_gradient(1,Val(2)), seed_gradient(1,Val(2)), false)
   x_log1p = MC{2,NS}(2.0, 3.0, Interval{Float64}(1.0,4.0), seed_gradient(1,Val(2)), seed_gradient(1,Val(2)), false)

   yref_log = MC{2,NS}(0.46209812037329695,1.09861, Interval{Float64}(0, 1.38629), SVector{2,Float64}([0.462098, 0.0]), SVector{2,Float64}([0.3333333, 0.0]), false)
   yref_log2 = MC{2,NS}(0.6666666666666666, 1.584962500721156, Interval{Float64}(0, 2), SVector{2,Float64}([0.666667, 0.0]), SVector{2,Float64}([0.4808983469629878, 0.0]), false)
   yref_log10 = MC{2,NS}(0.20068666377598746, 0.47712125471966244, Interval{Float64}(0, 0.60206), SVector{2,Float64}([0.200687, 0.0]), SVector{2,Float64}([0.14476482730108392, 0.0]), false)
   yref_log1p = MC{2,NS}(0.998577424517997, 1.3862943611198906, Interval{Float64}(0.693147, 1.60944), SVector{2,Float64}([0.30543, 0.0]), SVector{2,Float64}([0.25, 0.0]), false)

   @test check_vs_ref1(log, x_log, yref_log, mctol)
   @test check_vs_ref1(log2, x_log2, yref_log2, mctol)
   @test check_vs_ref1(log10, x_log10, yref_log10, mctol)
   @test check_vs_ref1(log1p, x_log1p, yref_log1p, mctol)

   #####  Square root #####
   x_sqrt_ns = MC{2,NS}(4.5, 4.5, Interval{Float64}(3.0,9.0), seed_gradient(1,Val(2)), seed_gradient(1,Val(2)), false)
   x_sqrt_d1 = MC{2,Diff}(4.0, 4.0, Interval{Float64}(3.0,7.0), seed_gradient(1,Val(2)), seed_gradient(1,Val(2)), false)

   yref_sqrt_ns = MC{2,NS}(2.049038105676658, 2.1213203435596424, Interval{Float64}(1.73205, 3.0), SVector{2,Float64}([0.211325, 0.0]), SVector{2,Float64}([0.235702, 0.0]), false)
   yref_sqrt_d1 = MC{2,Diff}(1.9604759334428057, 2.0, Interval{Float64}(1.73205, 2.64576), SVector{2,Float64}([0.228425, 0.0]), SVector{2,Float64}([0.25, 0.0]), false)

   @test check_vs_ref1(sqrt, x_sqrt_ns, yref_sqrt_ns, mctol)
   @test check_vs_ref1(sqrt, x_sqrt_d1, yref_sqrt_d1, mctol)

   #####  Step #####
   x_step_p_ns = MC{2,NS}(4.0, 4.0, Interval{Float64}(3.0,7.0), seed_gradient(1,Val(2)), seed_gradient(1,Val(2)), false)
   x_step_n_ns = MC{2,NS}(-4.0, -4.0, -Interval{Float64}(3.0,7.0), seed_gradient(1,Val(2)), seed_gradient(1,Val(2)), false)
   x_step_z_ns = MC{2,NS}(-2.0, -2.0, Interval{Float64}(-3.0,1.0), seed_gradient(1,Val(2)), seed_gradient(1,Val(2)), false)
   x_step_zp_ns = MC{2,NS}(0.5, 0.5, Interval{Float64}(-3.0,1.0), seed_gradient(1,Val(2)), seed_gradient(1,Val(2)), false)

   x_step_d1_p = MC{2,Diff}(4.0, 4.0, Interval{Float64}(3.0,7.0), seed_gradient(1,Val(2)), seed_gradient(1,Val(2)), false)
   x_step_d1_n = MC{2,Diff}(-4.0, -4.0, -Interval{Float64}(3.0,7.0), seed_gradient(1,Val(2)), seed_gradient(1,Val(2)), false)
   x_step_d1_z = MC{2,Diff}(-2.0, -2.0, Interval{Float64}(-3.0,1.0), seed_gradient(1,Val(2)), seed_gradient(1,Val(2)), false)

   yref_step_p = MC{2,NS}(1.0, 1.0, Interval{Float64}(1.0, 1.0), @SVector[0.0, 0.0], @SVector[0.0, 0.0], false)
   yref_step_n = MC{2,NS}(0.0, 0.0, Interval{Float64}(0.0, 0.0), @SVector[0.0, 0.0], @SVector[0.0, 0.0], false)
   yref_step_ns_z = MC{2,NS}(0.0, 0.3333333333333333, Interval{Float64}(0.0, 1.0), @SVector[0.0, 0.0], @SVector[-0.6666666666666666, 0.0], false)
   yref_step_ns_zp = MC{2,NS}(0.5, 1.0, Interval{Float64}(0.0, 1.0), @SVector[1.0, 0.0], @SVector[0.0, 0.0], false)

   yref_step_d1_p = MC{2,Diff}(1.0, 1.0, Interval{Float64}(1.0, 1.0), @SVector[0.0, 0.0], @SVector[0.0, 0.0], false)
   yref_step_d1_n = MC{2,Diff}(0.0, 0.0, Interval{Float64}(0.0, 0.0), @SVector[0.0, 0.0], @SVector[0.0, 0.0], false)
   yref_step_d1_z = MC{2,Diff}(0.0, 0.5555555555555556, Interval{Float64}(0.0, 1.0), @SVector[0.0, 0.0], @SVector[0.444444, 0.0], false)

   @test check_vs_ref1(step, x_step_p_ns, yref_step_p, mctol)
   @test check_vs_ref1(step, x_step_n_ns, yref_step_n, mctol)
   @test check_vs_ref1(step, x_step_z_ns, yref_step_ns_z, mctol)
   @test check_vs_ref1(step, x_step_zp_ns, yref_step_ns_zp, mctol)

   @test check_vs_ref1(step, x_step_d1_p, yref_step_d1_p, mctol)
   @test check_vs_ref1(step, x_step_d1_n, yref_step_d1_n, mctol)
   @test check_vs_ref1(step, x_step_d1_z, yref_step_d1_z, mctol)

   #####  Absolute value #####
   x_abs_ns = MC{2,NS}(4.5, 4.5, Interval{Float64}(-3.0,8.0), seed_gradient(1,Val(2)), seed_gradient(1,Val(2)), false)
   x_abs_d1 = MC{2,Diff}(2.0, 2.0, Interval{Float64}(-5.0,7.0), seed_gradient(1,Val(2)), seed_gradient(1,Val(2)), false)

   yref_abs_ns = MC{2,NS}(4.5, 6.409090909090908, Interval{Float64}(0.0, 8.0), @SVector[1.0, 0.0], @SVector[0.454545, 0.0], false)
   yref_abs_d1 = MC{2,Diff}(0.5714285714285714, 6.166666666666667, Interval{Float64}(0.0, 7.0), @SVector[0.5714285714285714, 0.0], @SVector[0.16666666666666666, 0.0], false)

   @test check_vs_ref1(abs, x_abs_ns, yref_abs_ns, mctol)
   @test check_vs_ref1(abs, x_abs_d1, yref_abs_d1, mctol)

   #####  Sign #####
   x_sign_d1_p = MC{2,Diff}(4.0, 4.0, Interval{Float64}(3.0,7.0), seed_gradient(1,Val(2)), seed_gradient(1,Val(2)), false)
   x_sign_d1_n = MC{2,Diff}(-4.0, -4.0, Interval{Float64}(-7.0,-3.0), seed_gradient(1,Val(2)), seed_gradient(1,Val(2)), false)
   x_sign_d1_z = MC{2,Diff}(-2.0, -2.0, Interval{Float64}(-3.0,1.0), seed_gradient(1,Val(2)), seed_gradient(1,Val(2)), false)

   yref_sign_d1_p = MC{2,Diff}(1.0, 1.0, Interval{Float64}(1.0, 1.0), @SVector[0.0, 0.0], @SVector[0.0, 0.0], false)
   yref_sign_d1_n = MC{2,Diff}(-1.0, -1.0, Interval{Float64}(-1.0, -1.0), @SVector[0.0, 0.0], @SVector[0.0, 0.0], false)
   yref_sign_d1_z = MC{2,Diff}(-1.0, 0.11111111111111116, Interval{Float64}(-1.0, 1.0), @SVector[0.0, 0.0], @SVector[0.888889, 0.0], false)

   @test check_vs_ref1(sign, x_sign_d1_p, yref_sign_d1_p, mctol)
   @test check_vs_ref1(sign, x_sign_d1_n, yref_sign_d1_n, mctol)
   @test check_vs_ref1(sign, x_sign_d1_z, yref_sign_d1_z, mctol)

   #####  Sine #####
   x_sin_p = MC{2,NS}(4.0, 4.0, Interval{Float64}(3.0,7.0), seed_gradient(1,Val(2)), seed_gradient(1,Val(2)), false)
   x_sin_n = MC{2,NS}(-4.0, -4.0, Interval{Float64}(-7.0,-3.0), seed_gradient(1,Val(2)), seed_gradient(1,Val(2)), false)
   x_sin_z = MC{2,NS}(-2.0, -2.0, Interval{Float64}(-3.0,1.0), seed_gradient(1,Val(2)), seed_gradient(1,Val(2)), false)
   x_sin_d1_z = MC{2,Diff}(-2.0, -2.0, Interval{Float64}(-3.0,1.0), seed_gradient(1,Val(2)), seed_gradient(1,Val(2)), false)

   yref_sin_ns_p = MC{2,NS}(-0.7568024953079283, 0.2700866557245978, Interval{Float64}(-1.0, 0.656987), @SVector[-0.653644, 0.0], @SVector[0.128967, 0.0], false)
   yref_sin_ns_n = MC{2,NS}(-0.2700866557245979, 0.7568024953079283, Interval{Float64}(-0.656987, 1.0), @SVector[0.128967, 0.0], @SVector[-0.653644, 0.0], false)
   yref_sin_ns_z = MC{2,NS}(-0.9092974268256817, 0.10452774015707458, Interval{Float64}(-1.0, 0.841471), @SVector[-0.416147, 0.0], @SVector[0.245648, 0.0], false)
   yref_sin_d1_z = MC{2,Diff}(-0.9092974268256817, 0.10452774015707458, Interval{Float64}(-1.0, 0.841471), @SVector[-0.416147, 0.0], @SVector[0.245648, 0.0], false)

   @test check_vs_ref1(sin, x_sin_p, yref_sin_ns_p, mctol)
   @test check_vs_ref1(sin, x_sin_n, yref_sin_ns_n, mctol)
   @test check_vs_ref1(sin, x_sin_z, yref_sin_ns_z, mctol)
   @test check_vs_ref1(sin, x_sin_d1_z, yref_sin_d1_z, mctol)

   #####  Cosine #####
   x_cos_p = MC{2,Diff}(4.0, 4.0, Interval{Float64}(3.0,7.0), seed_gradient(1,Val(2)), seed_gradient(1,Val(2)), false)
   x_cos_n = MC{2,Diff}(-4.0, -4.0, Interval{Float64}(-7.0,-3.0), seed_gradient(1,Val(2)), seed_gradient(1,Val(2)), false)
   x_cos_z = MC{2,Diff}(-2.0, -2.0, Interval{Float64}(-3.0,1.0), seed_gradient(1,Val(2)), seed_gradient(1,Val(2)), false)
   x_cos_z_ns = MC{2,NS}(-2.0, -2.0, Interval{Float64}(-3.0,1.0), seed_gradient(1,Val(2)), seed_gradient(1,Val(2)), false)

   yref_cos_d1_p = MC{2,Diff}(-0.703492113936536, -0.31034065427934965, Interval{Float64}(-1.0, 1.0), @SVector[0.485798, 0.0], @SVector[0.679652, 0.0], false)
   yref_cos_d1_n = MC{2,Diff}(-0.703492113936536, -0.31034065427934965, Interval{Float64}(-1.0, 1.0), @SVector[-0.485798, 0.0], @SVector[-0.679652, 0.0], false)
   yref_cos_d1_z = MC{2,Diff}(-0.6314158569813042, -0.222468094224762, Interval{Float64}(-0.989993, 1.0), @SVector[0.390573, 0.0], @SVector[0.76752, 0.0], false)
   yref_cos_ns_z = MC{2,NS}(-0.6314158569813042, -0.222468094224762, Interval{Float64}(-0.989993, 1.0), @SVector[0.390573, 0.0], @SVector[0.76752, 0.0], false)

   @test check_vs_ref1(cos, x_cos_p, yref_cos_d1_p, mctol)
   @test check_vs_ref1(cos, x_cos_n, yref_cos_d1_n, mctol)
   @test check_vs_ref1(cos, x_cos_z, yref_cos_d1_z, mctol)
   @test check_vs_ref1(cos, x_cos_z_ns, yref_cos_ns_z, mctol)

   ##### Tangent #####
   x_tan_p = MC{2,NS}(0.6, 0.6, Interval{Float64}(0.5,1.0), seed_gradient(1,Val(2)), seed_gradient(1,Val(2)), false)
   x_tan_n = MC{2,NS}(-0.8, -0.8, Interval{Float64}(-1.0,-0.5), seed_gradient(1,Val(2)), seed_gradient(1,Val(2)), false)
   x_tan_z = MC{2,NS}(-0.3,-0.3, Interval{Float64}(-0.5,0.5), seed_gradient(1,Val(2)), seed_gradient(1,Val(2)), false)
   x_tan_d1_z = MC{2,Diff}(-0.3,-0.3, Interval{Float64}(-0.5,0.5), seed_gradient(1,Val(2)), seed_gradient(1,Val(2)), false)
   x_tan_err = MC{2,Diff}(0.6, 0.6, Interval{Float64}(-4.5, 5.0), seed_gradient(1,Val(2)), seed_gradient(1,Val(2)), false)

   yref_tan_ns_p = MC{2,NS}(0.6841368083416923, 0.7485235368060128, Interval{Float64}(0.546302, 1.55741), @SVector[1.46804, 0.0], @SVector[2.02221, 0.0], false)
   yref_tan_ns_n = MC{2,NS}(-1.1529656307304577, -1.0296385570503641, Interval{Float64}(-1.55741, -0.546302), @SVector[2.02221, 0.0], @SVector[2.06016, 0.0], false)
   yref_tan_ns_z = MC{2,NS}(-0.332534, -0.30933624960962325, Interval{Float64}(-0.546303,0.546303), @SVector[1.06884, 0.0], @SVector[1.09569, 0.0], false)
   yref_tan_d1_z = MC{2,Diff}(-0.332534, -0.309336, Interval{Float64}(-0.546303,0.546303), @SVector[1.06884, 0.0], @SVector[1.09569, 0.0], false)

   @test check_vs_ref1(tan, x_tan_p, yref_tan_ns_p, mctol)
   @test check_vs_ref1(tan, x_tan_n, yref_tan_ns_n, mctol)
   @test check_vs_ref1(tan, x_tan_z, yref_tan_ns_z, mctol)
   @test check_vs_ref1(tan, x_tan_d1_z, yref_tan_d1_z, mctol)
   @test_throws ErrorException tan(x_tan_err)

   #####  Arcsine #####
   x_asin_p = MC{2,NS}(-0.7, -0.7, Interval{Float64}(-0.9,-0.5), seed_gradient(1,Val(2)), seed_gradient(1,Val(2)), false)
   x_asin_n = MC{2,NS}(0.7, 0.7, Interval{Float64}(0.5,0.9), seed_gradient(1,Val(2)), seed_gradient(1,Val(2)), false)
   x_asin_z = MC{2,NS}(-0.1, -0.1, Interval{Float64}(-0.5,0.5), seed_gradient(1,Val(2)), seed_gradient(1,Val(2)), false)
   x_asin_ns_z = MC{2,Diff}(-0.1, -0.1, Interval{Float64}(-0.5,0.5), seed_gradient(1,Val(2)), seed_gradient(1,Val(2)), false)

   yref_asin_d1_p = MC{2,NS}(-0.8216841452984665, -0.775397496610753, Interval{Float64}(-1.11977, -0.523598), @SVector[1.49043, 0.0], @SVector[1.40028, 0.0], false)
   yref_asin_d1_n = MC{2,NS}(0.775397496610753, 0.8216841452984665, Interval{Float64}(0.523598, 1.11977), @SVector[1.40028, 0.0], @SVector[1.49043, 0.0], false)
   yref_asin_d1_z = MC{2,NS}(-0.10958805193420748, -0.0974173098978382, Interval{Float64}(-0.523599, 0.523599), @SVector[1.03503, 0.0], @SVector[1.03503, 0.0], false)
   yref_asin_ns_z = MC{2,Diff}(-0.10958805193420748, -0.0974173098978382, Interval{Float64}(-0.523599, 0.523599), @SVector[1.03503, 0.0], @SVector[1.03503, 0.0], false)

   @test check_vs_ref1(asin, x_asin_p, yref_asin_d1_p, mctol)
   @test check_vs_ref1(asin, x_asin_n, yref_asin_d1_n, mctol)
   @test check_vs_ref1(asin, x_asin_z, yref_asin_d1_z, mctol)
   @test check_vs_ref1(asin, x_asin_ns_z, yref_asin_ns_z, mctol)

   ##### Arctangent #####
   x_atan_p = MC{2,NS}(4.0, 4.0, Interval{Float64}(3.0, 7.0), seed_gradient(1,Val(2)), seed_gradient(1,Val(2)), false)
   yref_atan_d1_p = MC{2,NS}(1.294009147346374, 1.3258176636680326, Interval{Float64}(1.24904, 1.4289), @SVector[0.04496337494811958, 0.0], @SVector[0.058823529411764705, 0.0], false)
   @test check_vs_ref1(atan, x_atan_p, yref_atan_d1_p, mctol)

   ##### Hyperbolic Sine #####
   x_sinh_p = MC{2,Diff}(4.0, 4.0, Interval{Float64}(3.0,7.0), seed_gradient(1,Val(2)), seed_gradient(1,Val(2)), false)
   x_sinh_n = MC{2,Diff}(-4.0, -4.0, Interval{Float64}(-7.0,-3.0), seed_gradient(1,Val(2)), seed_gradient(1,Val(2)), false)
   x_sinh_z = MC{2,Diff}(-2.0, -2.0, Interval{Float64}(-3.0,1.0), seed_gradient(1,Val(2)), seed_gradient(1,Val(2)), false)
   x_sinh_ns_p = MC{2,NS}(4.0, 4.0, Interval{Float64}(3.0,7.0), seed_gradient(1,Val(2)), seed_gradient(1,Val(2)), false)

   yref_sinh_d1_p = MC{2,Diff}(27.28991719712775, 144.59243701386904, Interval{Float64}(10.0178, 548.3161232732466), @SVector[27.3082, 0.0], @SVector[134.57456208645914, 0.0], false)
   yref_sinh_d1_n = MC{2,Diff}(-144.59243701386904, -27.28991719712775, Interval{Float64}(-548.3161232732466, -10.0178), @SVector[134.57456208645914 , 0.0], @SVector[27.3082, 0.0], false)
   yref_sinh_d1_z = MC{2,Diff}(-6.212568527712605, -3.626860407847019, Interval{Float64}(-10.0179, 1.17521), @SVector[3.8053063996972973, 0.0], @SVector[3.7622, 0.0], false)
   yref_sinh_ns_p = MC{2,NS}(27.28991719712775, 144.59243701386904, Interval{Float64}(10.0178, 548.317), @SVector[27.3082, 0.0], @SVector[134.575, 0.0], false)

   @test check_vs_ref1(sinh, x_sinh_p, yref_sinh_d1_p, mctol)
   @test check_vs_ref1(sinh, x_sinh_n, yref_sinh_d1_n, mctol)
   @test check_vs_ref1(sinh, x_sinh_z, yref_sinh_d1_z, mctol)
   @test_broken check_vs_ref1(sinh, x_sinh_ns_p, yref_sinh_ns_p, mctol)

   ##### Hyperbolic Cosine #####
   x_cosh_ns = MC{2,NS}(4.0, 4.0, Interval{Float64}(3.0,7.0), seed_gradient(1,Val(2)), seed_gradient(1,Val(2)), false)
   x_cosh_d1 = MC{2,Diff}(4.0, 4.0, Interval{Float64}(3.0,7.0), seed_gradient(1,Val(2)), seed_gradient(1,Val(2)), false)
   yref_cosh_ns = MC{2,NS}(27.308232836016487, 144.63000528563632, Interval{Float64}(10.0676, 548.317035), @SVector[27.28991719712775, 0.0], @SVector[134.56234328985857, 0.0], false)
   yref_cosh_d1 = MC{2,Diff}(27.308232836016487, 144.63000528563632, Interval{Float64}(10.0676, 548.317035), @SVector[27.28991719712775, 0.0], @SVector[134.56234328985857, 0.0], false)
   @test check_vs_ref1(cosh, x_cosh_ns, yref_cosh_ns, mctol)
   @test check_vs_ref1(cosh, x_cosh_d1, yref_cosh_d1, mctol)

   ##### Hyperbolic Tangent #####
   x_tanh_p = MC{2,Diff}(4.0, 4.0, Interval{Float64}(3.0,7.0), seed_gradient(1,Val(2)), seed_gradient(1,Val(2)), false)
   x_tanh_n = MC{2,Diff}(-4.0, -4.0, Interval{Float64}(-7.0,-3.0), seed_gradient(1,Val(2)), seed_gradient(1,Val(2)), false)
   x_tanh_z1 = MC{2,Diff}(-2.0, -2.0, Interval{Float64}(-3.0,1.0), seed_gradient(1,Val(2)), seed_gradient(1,Val(2)), false)
   x_tanh_z2 = MC{2,Diff}(2.0, 2.0, Interval{Float64}(-1.0,3.0), seed_gradient(1,Val(2)), seed_gradient(1,Val(2)), false)
   x_tanh_z1_ns = MC{2,NS}(-2.0, -2.0, Interval{Float64}(-3.0,1.0), seed_gradient(1,Val(2)), seed_gradient(1,Val(2)), false)

   yref_tanh_d1_p = MC{2,Diff}(0.996290649501034, 0.999329299739067, Interval{Float64}(0.995054,0.999999), @SVector[0.0012359, 0.0], @SVector[0.00134095, 0.0], false)
   yref_tanh_d1_n = MC{2,Diff}(-0.999329299739067, -0.996290649501034, Interval{Float64}(-0.999999, -0.995054), @SVector[0.00134095, 0.0], @SVector[0.0012359, 0.0], false)
   yref_tanh_d1_z1 = MC{2,Diff}(-0.9640275800758169, -0.5558207301372651, Interval{Float64}(-0.995055, 0.761595), @SVector[0.0706508, 0.0], @SVector[0.439234, 0.0], false)
   yref_tanh_d1_z2 = MC{2,Diff}(0.5558207301372651, 0.9640275800758169, Interval{Float64}(-0.761595, 0.995055), @SVector[0.439234, 0.0], @SVector[0.0706508, 0.0], false)
   yref_tanh_ns = MC{2,NS}(-0.9640275800758169, -0.5558207301372651, Interval{Float64}(-0.995055, 0.761595), @SVector[0.0706508, 0.0], @SVector[0.439234, 0.0], false)

   @test check_vs_ref1(tanh, x_tanh_p, yref_tanh_d1_p, mctol)
   @test check_vs_ref1(tanh, x_tanh_n, yref_tanh_d1_n, mctol)
   @test check_vs_ref1(tanh, x_tanh_z1, yref_tanh_d1_z1, mctol)
   @test check_vs_ref1(tanh, x_tanh_z2, yref_tanh_d1_z2, mctol)
   @test check_vs_ref1(tanh, x_tanh_z1_ns, yref_tanh_ns, mctol)

   ##### Inverse Hyperbolic Sine #####
   x_asinh_p = MC{2,Diff}(0.3, 0.3, Interval{Float64}(0.1,0.7), seed_gradient(1,Val(2)), seed_gradient(1,Val(2)), false)
   x_asinh_n = MC{2,Diff}(-0.3, -0.3, Interval{Float64}(-0.7,-0.1), seed_gradient(1,Val(2)), seed_gradient(1,Val(2)), false)
   x_asinh_z1 = MC{2,Diff}(2.0, 2.0, Interval{Float64}(-3.0,7.0), seed_gradient(1,Val(2)), seed_gradient(1,Val(2)), false)
   x_asinh_z2 = MC{2,Diff}(-2.0, -2.0, Interval{Float64}(-7.0,3.0), seed_gradient(1,Val(2)), seed_gradient(1,Val(2)), false)
   x_asinh_z1_ns = MC{2,NS}(2.0, 2.0, Interval{Float64}(-3.0,7.0), seed_gradient(1,Val(2)), seed_gradient(1,Val(2)), false)

   yref_asinh_d1_p = MC{2,Diff}(0.2841115746269236, 0.29567304756342244, Interval{Float64}(0.099834,0.652667), @SVector[0.921387, 0.0], @SVector[0.957826, 0.0], false)
   yref_asinh_d1_n = MC{2,Diff}(-0.29567304756342244, -0.2841115746269236, Interval{Float64}(-0.652667,-0.099834), @SVector[0.957826, 0.0], @SVector[0.921387, 0.0], false)
   yref_asinh_d1_z1 = MC{2,Diff}(0.3730697449603356, 1.4436354751788103, Interval{Float64}(-1.81845, 2.64413), @SVector[0.45421, 0.0], @SVector[0.447214, 0.0], false)
   yref_asinh_d1_z2 = MC{2,Diff}(-1.4436354751788103, -0.3730697449603356, Interval{Float64}(-2.64413,1.81845), @SVector[0.447214, 0.0], @SVector[0.45421, 0.0], false)
   yref_asinh_ns = MC{2,NS}(0.37306974496033596, 1.4436354751788103, Interval{Float64}(-1.8184464592320668, 2.6441207610586295), @SVector[0.45421020321965866, 0.0], @SVector[0.44721359549, 0.0], false)

   @test check_vs_ref1(asinh, x_asinh_p, yref_asinh_d1_p, mctol)
   @test check_vs_ref1(asinh, x_asinh_n, yref_asinh_d1_n, mctol)
   @test check_vs_ref1(asinh, x_asinh_z1, yref_asinh_d1_z1, mctol)
   @test check_vs_ref1(asinh, x_asinh_z2, yref_asinh_d1_z2, mctol)
   @test check_vs_ref1(asinh, x_asinh_z1_ns, yref_asinh_ns, mctol)

   ##### Inverse Hyperbolic Cosine #####
   x_acosh = MC{2,NS}(4.0, 4.0, Interval{Float64}(3.0,7.0), seed_gradient(1,Val(2)), seed_gradient(1,Val(2)), false)
   yref_acosh = MC{2,NS}(1.9805393289917226, 2.0634370688955608, Interval{Float64}(1.76274,2.63392), @SVector[0.217792, 0.0], @SVector[0.258199, 0.0], false)
   @test check_vs_ref1(acosh, x_acosh, yref_acosh, mctol)

   ##### Inverse Hyperbolic Tangent #####
   x_atanh_p = MC{2,Diff}(0.6, 0.6, Interval{Float64}(0.1,0.7), seed_gradient(1,Val(2)), seed_gradient(1,Val(2)), false)
   x_atanh_n = MC{2,Diff}(-0.6, -0.6, Interval{Float64}(-0.7,-0.1), seed_gradient(1,Val(2)), seed_gradient(1,Val(2)), false)
   x_atanh_z1 = MC{2,Diff}(0.6, 0.6, Interval{Float64}(-0.6,0.7), seed_gradient(1,Val(2)), seed_gradient(1,Val(2)), false)
   x_atanh_z2 = MC{2,Diff}(-0.5, -0.5, Interval{Float64}(-0.7,0.6), seed_gradient(1,Val(2)), seed_gradient(1,Val(2)), false)
   x_atanh_z1_ns = MC{2,NS}(0.6, 0.6, Interval{Float64}(-0.6,0.7), seed_gradient(1,Val(2)), seed_gradient(1,Val(2)), false)

   yref_atanh_d1_p = MC{2,Diff}(0.6931471805599453, 0.7394729977002236, Interval{Float64}(0.100335,0.867301), @SVector[1.5625, 0.0], @SVector[1.27828, 0.0], false)
   yref_atanh_d1_n = MC{2,Diff}(-0.7394729977002236, -0.6931471805599453, Interval{Float64}(-0.867301,-0.100335), @SVector[1.27828, 0.0], @SVector[1.5625, 0.0], false)
   yref_atanh_d1_z1 = MC{2,Diff}(0.6931471805599453, 0.7263562565915014, Interval{Float64}(-0.6931471805599453,0.867301), @SVector[1.5625, 0.0], @SVector[1.4094427110255183, 0.0], false)
   yref_atanh_d1_z2 = MC{2,Diff}(-0.5854119854889636, -0.5493061443340549, Interval{Float64}(-0.867301,0.6931471805599453), @SVector[1.4094427110254484, 0.0], @SVector[1.3333333333333333, 0.0], false)
   yref_atanh_ns = MC{2,NS}(0.6931471805599453, 0.72635625, Interval{Float64}(-0.6931471805599453,0.867301), @SVector[1.5625, 0.0], @SVector[1.4094427110255183, 0.0], false)

   @test check_vs_ref1(atanh, x_atanh_p, yref_atanh_d1_p, mctol)
   @test check_vs_ref1(atanh, x_atanh_n, yref_atanh_d1_n, mctol)
   @test check_vs_ref1(atanh, x_atanh_z1, yref_atanh_d1_z1, mctol)
   @test check_vs_ref1(atanh, x_atanh_z2, yref_atanh_d1_z2, mctol)
   @test check_vs_ref1(atanh, x_atanh_z1_ns, yref_atanh_ns, mctol)

   # CONVERSION
   X = MC{2,NS}(4.5, 4.5, Interval{Float64}(-3.0,8.0), seed_gradient(1,Val(2)), seed_gradient(1,Val(2)), false)
   X1 = convert(MC{2,NS}, 1)
   X2 = convert(MC{2,NS}, 1.1)
   X3 = convert(MC{2,NS}, Interval{Float64}(2.1,4.3))
   @test X1.cc == 1.0
   @test X1.cv == 1.0
   @test X2.cc == 1.1
   @test X2.cv == 1.1
   @test X3.cc == 4.3
   @test X3.cv == 2.1

   @test +X == X
end


@testset "Test Multivariant w/Constant" begin

   mctol = 1E-4

   x = MC{2,NS}(4.5, 4.5, Interval{Float64}(-3.0,8.0), seed_gradient(1,Val(2)), seed_gradient(1,Val(2)), false)

   x1 = x + 2.1; @test x1.cv == x1.cc == 6.6
   x2 = 2.3 + x; @test x2.cv == x2.cc == 6.8
   x3 = x + 2;   @test x3.cv == x3.cc == 6.5
   x4 = 3 + x;   @test x4.cv == x4.cc == 7.5

   x1 = x + Float16(2.1); @test x1.cv == x1.cc == 6.599609375
   x2 = Float16(2.3) + x; @test x2.cv == x2.cc == 6.80078125
   x3 = x + Int16(2);     @test x3.cv == x3.cc == 6.5
   x4 = Int16(3) + x;     @test x4.cv == x4.cc == 7.5

   x1 = x - 2.1; @test x1.cv == x1.cc == 2.4
   x2 = 2.3 - x; @test x2.cv == x2.cc == -2.2
   x3 = x - 2;   @test x3.cv == x3.cc == 2.5
   x4 = 3 - x;   @test x4.cv == x4.cc == -1.5

   x1 = x - Float16(2.1); @test x1.cv == x1.cc == 2.400390625
   x2 = Float16(2.3) - x; @test x2.cv == x2.cc == -2.19921875
   x3 = x - Int16(2);     @test x3.cv == x3.cc == 2.5
   x4 = Int16(3) - x;     @test x4.cv == x4.cc == -1.5

   x1 = x * 2.1; @test x1.cv == x1.cc == 9.450000000000001
   x2 = 2.3 * x; @test x2.cv == x2.cc == 10.35
   x3 = x * 2;   @test x3.cv == x3.cc == 9.0
   x4 = 3 * x;   @test x4.cv == x4.cc == 13.5

   x1 = x * Float16(2.1); @test x1.cv == x1.cc == 9.4482421875
   x2 = Float16(2.3) * x; @test x2.cv == x2.cc == 10.353515625
   x3 = x * Int16(2);     @test x3.cv == x3.cc == 9.0
   x4 =  Int16(3) * x;    @test x4.cv == x4.cc == 13.5

   x1 = x * (-2.1); @test x1.cv == x1.cc == -9.450000000000001
   x2 = (-2.3) * x; @test x2.cv == x2.cc == -10.35
   x3 = x * (-2);   @test x3.cv == x3.cc == -9.0
   x4 = (-3) * x;   @test x4.cv == x4.cc == -13.5

   x1 = x * Float16(-2.1); @test x1.cv == x1.cc == -9.4482421875
   x2 = Float16(-2.3) * x; @test x2.cv == x2.cc == -10.353515625
   x3 = x * Int16(-2);     @test x3.cv == x3.cc == -9.0
   x4 =  Int16(-3) * x;    @test x4.cv == x4.cc == -13.5

   x1 = 1.0/(x+4)
   @test isapprox(x1.cc, 0.37500000000000017, atol=mctol)
   @test isapprox(x1.cv, 0.11764705882352941, atol=mctol)
   @test isapprox(x1.cc_grad[1], -0.08333333333333334, atol=mctol)
   @test isapprox(x1.cc_grad[2], -0.0, atol=mctol)
   @test isapprox(x1.cv_grad[1], -0.01384083044982699, atol=mctol)
   @test isapprox(x1.cv_grad[2], -0.0, atol=mctol)

   x2 = 1/(x+4)
   @test isapprox(x2.cc, 0.37500000000000017, atol=1E-6)
   @test isapprox(x2.cv, 0.11764705882352941, atol=1E-6)
   @test isapprox(x2.cc_grad[1], -0.08333333333333334, atol=1E-6)
   @test isapprox(x2.cc_grad[2], -0.0, atol=1E-6)
   @test isapprox(x2.cv_grad[1], -0.01384083044982699, atol=1E-6)
   @test isapprox(x2.cv_grad[2], -0.0, atol=1E-6)

   x_min_ns = MC{2,NS}(1.0909090909090908, 3.0, Interval{Float64}(-3.0, 3.0), @SVector[0.545455, 0.0], @SVector[0.0, 0.0], false)
   x_max_ns = MC{2,NS}(1.0909090909090908, 3.0, Interval{Float64}(-3.0, 3.0), @SVector[0.545455, 0.0], @SVector[0.0, 0.0], false)
   x_min_d1 = MC{2,Diff}(1.0909090909090908, 3.0, Interval{Float64}(-3.0, 3.0), @SVector[0.545455, 0.0], @SVector[0.0, 0.0], false)
   x_max_d1 = MC{2,Diff}(1.0909090909090908, 3.0, Interval{Float64}(-3.0, 3.0), @SVector[0.545455, 0.0], @SVector[0.0, 0.0], false)

   yref1_min_ns = MC{2,NS}(1.0909090909090908, 3.0, Interval{Float64}(-3.0, 3.0), @SVector[0.545455, 0.0], @SVector[0.0, 0.0], false)
   yref2_min_ns = MC{2,NS}(1.0909090909090908, 3.0, Interval{Float64}(-3.0, 3.0), @SVector[0.545455, 0.0], @SVector[0.0, 0.0], false)
   yref1_max_ns = MC{2,NS}(3.0, 3.0, Interval{Float64}(3.0, 3.0), @SVector[0.0, 0.0], @SVector[0.0, 0.0], false)
   yref2_max_ns = MC{2,NS}(3.0, 3.0, Interval{Float64}(3.0, 3.0), @SVector[0.0, 0.0], @SVector[0.0, 0.0], false)
   yref1_min_d1 = MC{2,Diff}(1.0909090909090908, 3.0, Interval{Float64}(-3.0, 3.0), @SVector[0.545455, 0.0], @SVector[0.0, 0.0], false)
   yref2_min_d1 = MC{2,Diff}(1.0909090909090908, 3.0, Interval{Float64}(-3.0, 3.0), @SVector[0.545455, 0.0], @SVector[0.0, 0.0], false)
   yref1_max_d1 = MC{2,Diff}(3.0, 3.0, Interval{Float64}(3.0, 3.0), @SVector[0.0, 0.0], @SVector[0.0, 0.0], false)
   yref2_max_d1 = MC{2,Diff}(3.0, 3.0, Interval{Float64}(3.0, 3.0), @SVector[0.0, 0.0], @SVector[0.0, 0.0], false)

   @test check_vs_refv(min, x_min_ns, 3.0, yref1_min_ns, mctol)
   @test check_vs_refv(min, x_min_ns, 3.0, yref2_min_ns, mctol)
   @test check_vs_refv(max, x_max_ns, 3.0, yref1_max_ns, mctol)
   @test check_vs_refv(max, x_max_ns, 3.0, yref2_max_ns, mctol)

   @test check_vs_refv(min, x_min_d1, 3.0, yref1_min_d1, mctol)
   @test check_vs_refv(min, x_min_d1, 3.0, yref2_min_d1, mctol)
   @test check_vs_refv(max, x_max_d1, 3.0, yref1_max_d1, mctol)
   @test check_vs_refv(max, x_max_d1, 3.0, yref2_max_d1, mctol)

   x1 = MC{2,NS}(4.0, 4.0, Interval{Float64}(3.0, 7.0), seed_gradient(1,Val(2)), seed_gradient(1,Val(2)), false)
   x2 = MC{2,NS}(-4.5, -4.5, Interval{Float64}(-8.0, -3.0), seed_gradient(1,Val(2)), seed_gradient(1,Val(2)), false)
   x3 = MC{2,NS}(4.5, 4.5, Interval{Float64}(3.0, 8.0), seed_gradient(1,Val(2)), seed_gradient(1,Val(2)), false)
   x4 = MC{2,NS}(-4.5, -4.5, Interval{Float64}(-8.0, -3.0), seed_gradient(1,Val(2)), seed_gradient(1,Val(2)), false)
   x5 = MC{2,NS}(-4.0, -4.0, Interval{Float64}(-7.0, -3.0), seed_gradient(1,Val(2)), seed_gradient(1,Val(2)), false)
   x6 = MC{2,NS}(-2.0, -2.0,Interval{Float64}(-3.0, 1.0), seed_gradient(1,Val(2)), seed_gradient(1,Val(2)), false)

   yref_ns_pow1 = MC{2,NS}(16.0, 19.0, Interval{Float64}(9.0, 49.0), @SVector[8.0, 0.0], @SVector[10.0, 0.0], false)
   yref_ns_pow2 = MC{2,NS}(0.0625, 0.08843537414965986, Interval{Float64}(0.0204081, 0.111112), @SVector[-0.03125, 0.0], @SVector[-0.0226757, 0.0], false)
   yref_ns_pow3 = MC{2,NS}(0.04938271604938271, 0.08246527777777776, Interval{Float64}(0.015625, 0.111112), @SVector[0.0219479, 0.0], @SVector[0.0190972, 0.0], false)
   yref_ns_pow4 = MC{2,NS}(20.25, 25.5, Interval{Float64}(9.0, 64.0), @SVector[-9.0, 0.0], @SVector[-11.0, 0.0], false)
   yref_ns_pow5 = MC{2,NS}(-172.5, -91.125, Interval{Float64}(-512.0, -27.0), @SVector[97.0, 0.0], @SVector[60.75, 0.0], false)
   yref_ns_pow6 = MC{2,NS}(410.0625, 1285.5, Interval{Float64}(81.0, 4096.0), @SVector[-364.5, 0.0], @SVector[-803.0, 0.0], false)
   yref_ns_pow7 = MC{2,NS}(91.125, 172.5, Interval{Float64}(27.0, 512.0), @SVector[60.75, 0.0], @SVector[97.0, 0.0], false)
   yref_ns_pow8 = MC{2,NS}(410.0625, 1285.5, Interval{Float64}(81.0, 4096.0), @SVector[364.5, 0.0], @SVector[803.0, 0.0], false)
   yref_ns_pow9 = MC{2,NS}(410.0625, 1285.5, Interval{Float64}(81.0, 4096.0), @SVector[-364.5, 0.0], @SVector[-803.0, 0.0], false)
   @test check_vs_refv(^, x1, 2, yref_ns_pow1, mctol)
   @test check_vs_refv(^, x1, -2, yref_ns_pow2, mctol)
   @test check_vs_refv(^, x2, -2, yref_ns_pow3, mctol)
   @test check_vs_refv(^, x2, 2, yref_ns_pow4, mctol)
   @test check_vs_refv(^, x2, 3, yref_ns_pow5, mctol)
   @test check_vs_refv(^, x2, 4, yref_ns_pow6, mctol)
   @test check_vs_refv(^, x3, 3, yref_ns_pow7, mctol)
   @test check_vs_refv(^, x3, 4, yref_ns_pow8, mctol)
   @test check_vs_refv(^, x4, 4, yref_ns_pow9, mctol)
   @test x2^(1) == x2

   x1_d1 = MC{2,Diff}(4.0, 4.0, Interval{Float64}(3.0, 7.0), seed_gradient(1,Val(2)), seed_gradient(1,Val(2)), false)
   x2_d1 = MC{2,Diff}(-4.5, -4.5, Interval{Float64}(-8.0, -3.0), seed_gradient(1,Val(2)), seed_gradient(1,Val(2)), false)
   x3_d1 = MC{2,Diff}(4.5, 4.5, Interval{Float64}(3.0, 8.0), seed_gradient(1,Val(2)), seed_gradient(1,Val(2)), false)
   x4_d1 = MC{2,Diff}(-4.5, -4.5, Interval{Float64}(-8.0, -3.0), seed_gradient(1,Val(2)), seed_gradient(1,Val(2)), false)
   x5_d1 = MC{2,Diff}(-4.0, -4.0, Interval{Float64}(-7.0, -3.0), seed_gradient(1,Val(2)), seed_gradient(1,Val(2)), false)
   x6_d1 = MC{2,Diff}(-2.0, -2.0,Interval{Float64}(-3.0, 1.0), seed_gradient(1,Val(2)), seed_gradient(1,Val(2)), false)

   yref_d1_pow1 = MC{2,Diff}(0.25, 0.2857142857142857, Interval{Float64}(0.142857, 0.333334), @SVector[-0.0625, 0.0], @SVector[-0.047619, 0.0], false)
   yref_d1_pow2 = MC{2,Diff}(16.0, 19.0, Interval{Float64}(9.0, 49.0), @SVector[8.0, 0.0], @SVector[10.0, 0.0], false)
   yref_d1_pow3 = MC{2,Diff}(16.0, 19.0, Interval{Float64}(9.0, 49.0), @SVector[-8.0, 0.0], @SVector[-10.0, 0.0], false)
   yref_d1_pow4 = MC{2,Diff}(2.66666666666666, 7.0, Interval{Float64}(0.0, 9.0), @SVector[-4.0, 0.0], @SVector[-2.0, 0.0], false)
   yref_d1_pow5 = MC{2,Diff}(64.0, 106.0, Interval{Float64}(27.0, 343.0), @SVector[48.0, 0.0], @SVector[79.0, 0.0], false)
   yref_d1_pow6 = MC{2,Diff}(-106.0, -64.0, Interval{Float64}(-343.0, -27.0), @SVector[79.0, 0.0], @SVector[48.0, 0.0], false)
   yref_d1_pow7 = MC{2,Diff}(-20.25, -7.750, Interval{Float64}(-27.0, 1.0), @SVector[6.75, 0.0], @SVector[12.25, 0.0], false)
   yref_d1_pow8 = MC{2,Diff}(0.015625, 0.02850664075153871, Interval{Float64}(0.00291545, 0.0370371), @SVector[-0.0117188, 0.0], @SVector[-0.0085304, 0.0], false)
   yref_d1_pow9 = MC{2,Diff}(-0.02850664075153871, -0.015625, Interval{Float64}(-0.0370371, -0.00291545), @SVector[-0.0085304, 0.0], @SVector[-0.0117188, 0.0], false)
   yref_d1_pow10 = MC{2,Diff}(0.00390625, 0.009363382541225106, Interval{Float64}(0.000416493, 0.0123457), @SVector[-0.00390625, 0.0], @SVector[-0.0029823, 0.0], false)
   yref_d1_pow11 = MC{2,Diff}(0.00390625, 0.009363382541225106, Interval{Float64}(0.000416493, 0.0123457), @SVector[0.00390625, 0.0], @SVector[0.0029823, 0.0], false)
   yref_d1_pow12 = MC{2,Diff}(256.0, 661.0, Interval{Float64}(81.0, 2401.0), @SVector[256.0, 0.0], @SVector[580.0, 0.0], false)
   yref_d1_pow13 = MC{2,Diff}(256.0, 661.0, Interval{Float64}(81.0, 2401.0), @SVector[-256.0, 0.0], @SVector[-580.0, 0.0], false)
   yref_d1_pow14 = MC{2,Diff}(16.0, 61.0,  Interval{Float64}(0.0, 81.0), @SVector[-32.0, 0.0], @SVector[-20.0, 0.0], false)
   @test check_vs_refv(^, x1_d1, -1, yref_d1_pow1, mctol)
   @test check_vs_refv(^, x1_d1, 2, yref_d1_pow2, mctol)
   @test check_vs_refv(^, x5_d1, 2, yref_d1_pow3, mctol)
   @test check_vs_refv(^, x6_d1, 2, yref_d1_pow4, mctol)
   @test check_vs_refv(^, x1_d1, 3, yref_d1_pow5, mctol)
   @test check_vs_refv(^, x5_d1, 3, yref_d1_pow6, mctol)
   @test check_vs_refv(^, x6_d1, 3, yref_d1_pow7, mctol)
   @test check_vs_refv(^, x1_d1, -3, yref_d1_pow8, mctol)
   @test check_vs_refv(^, x5_d1, -3, yref_d1_pow9, mctol)
   @test check_vs_refv(^, x1_d1, -4, yref_d1_pow10, mctol)
   @test check_vs_refv(^, x5_d1, -4, yref_d1_pow11, mctol)
   @test check_vs_refv(^, x1_d1, 4, yref_d1_pow12, mctol)
   @test check_vs_refv(^, x5_d1, 4, yref_d1_pow13, mctol)
   @test check_vs_refv(^, x6_d1, 4, yref_d1_pow14, mctol)
end

@testset "Multiplication Operator" begin

    mctol = 1E-4

    ##### Case 1 #####
    x1 = MC{2,NS}(0.0, 0.0, Interval{Float64}(-2.0,1.0), seed_gradient(1,Val(2)), seed_gradient(1,Val(2)), false)
    y1 = MC{2,NS}(1.0, 1.0, Interval{Float64}(-1.0,2.0), seed_gradient(2,Val(2)), seed_gradient(2,Val(2)), false)
    yref1 = MC{2,NS}(-1.0, 2.0, Interval{Float64}(-4.0,2.0), @SVector[2.0, 1.0], @SVector[2.0, -2.0], false)
    @test check_vs_ref2(*, x1, y1, yref1, mctol)

    ##### Case 2 #####
    x2 = MC{2,NS}(3.0, 3.0, Interval{Float64}(1.0,5.0), seed_gradient(1,Val(2)), seed_gradient(1,Val(2)), false)
    y2 = MC{2,NS}(1.0, 1.0, Interval{Float64}(-1.0,2.0), seed_gradient(2,Val(2)), seed_gradient(2,Val(2)), false)
    yref2 = MC{2,NS}(1.0, 5.0, Interval{Float64}(-5.0,10.0), @SVector[2.0, 5.0], @SVector[2.0, 1.0], false)
    @test check_vs_ref2(*, x2, y2, yref2, mctol)

    ##### Case 3 #####
    x3 = MC{2,NS}(-4.0, -4.0, Interval{Float64}(-6.0,-2.0), seed_gradient(1,Val(2)), seed_gradient(1,Val(2)), false)
    y3 = MC{2,NS}(2.0, 2.0, Interval{Float64}(1.0,3.0), seed_gradient(2,Val(2)), seed_gradient(2,Val(2)), false)
    yref3 = MC{2,NS}(-10.0, -6.0, Interval{Float64}(-18.000000000000004,-1.9999999999999998), @SVector[3.0, -2.0], @SVector[1.0, -2.0], false)
    @test check_vs_ref2(*, x3, y3, yref3, mctol)

    ##### Case 4 #####
    x4 = MC{2,NS}(-4.0, -4.0, Interval{Float64}(-6.0,-2.0), seed_gradient(1,Val(2)), seed_gradient(1,Val(2)), false)
    y4 = MC{2,NS}(-5.0, -5.0, Interval{Float64}(-7.0,-3.0), seed_gradient(2,Val(2)), seed_gradient(2,Val(2)), false)
    yref4 = MC{2,NS}(16.0, 24.0, Interval{Float64}(6.0,42.0), @SVector[-3.0, -2.0], @SVector[-7.0, -2.0], false)
    @test check_vs_ref2(*, x4, y4, yref4, mctol)

    ##### Case 5 #####
    x5 = MC{2,NS}(-4.0, -4.0, Interval{Float64}(-6.0,-2.0), seed_gradient(1,Val(2)), seed_gradient(1,Val(2)), false)
    y5 = MC{2,NS}(-5.0, -5.0, Interval{Float64}(-7.0,4.0), seed_gradient(2,Val(2)), seed_gradient(2,Val(2)), false)
    yref5 = MC{2,NS}(16.0, 24.0, Interval{Float64}(-24.000000000000004, 42.00000000000001), @SVector[-7.0, -6.0], @SVector[-7.0, -2.0], false)
    @test check_vs_ref2(*, x5, y5, yref5, mctol)

    ##### Case 6 #####
    x6 = MC{2,NS}(-2.0, -2.0, Interval{Float64}(-3.0,4.0), seed_gradient(1,Val(2)), seed_gradient(1,Val(2)), false)
    y6 = MC{2,NS}(3.0, 3.0, Interval{Float64}(1.0,4.0), seed_gradient(2,Val(2)), seed_gradient(2,Val(2)), false)
    yref6 = MC{2,NS}(-8.0, -5.0, Interval{Float64}(-12.0, 16.0), @SVector[1.0, -3.0], @SVector[4.0, -3.0], false)
    @test check_vs_ref2(*, x6, y6, yref6, mctol)

    ##### Case 7 #####
    x7 = MC{2,NS}(-2.0, -2.0, Interval{Float64}(-3.0,4.0), seed_gradient(1,Val(2)), seed_gradient(1,Val(2)), false)
    y7 = MC{2,NS}(-4.0, -4.0, Interval{Float64}(-5.0,-3.0), seed_gradient(2,Val(2)), seed_gradient(2,Val(2)), false)
    yref7 = MC{2,NS}(7.0, 9.0, Interval{Float64}(-20.0, 15.0), @SVector[-5.0, -3.0], @SVector[-3.0, -3.0], false)
    @test check_vs_ref2(*, x7, y7, yref7, mctol)

    ##### Testing for Smooth Standard Mult #####
    seed1 = seed_gradient(1,Val(2))
    seed2 = seed_gradient(2,Val(2))
    x1 = MC{2,Diff}(0.0,0.0,Interval{Float64}(-200.0,200.0),seed1,seed1,false)
    y1 = MC{2,Diff}(200.0,200.0,Interval{Float64}(0.0,400.0),seed2,seed2,false)
    z1 = x1*y1
    @test isapprox(z1.cc,40000,atol=1E-4)
    @test isapprox(z1.cv,-40000,atol=1E-4)

    x2 = MC{2,Diff}(170.0,170.0,Interval{Float64}(100.0,240.0),seed1,seed1,false)
    y2 = MC{2,Diff}(250.0,250.0,Interval{Float64}(100.0,400.0),seed2,seed2,false)
    z2 = x2*y2
    @test isapprox(z2.cc,53000,atol=1E-4)
    @test isapprox(z2.cv,32000,atol=1E-4)

    x3 = MC{2,Diff}(-200.0,-200.0,Interval{Float64}(-300.0,-100.0),seed1,seed1,false)
    y3 = MC{2,Diff}(-300.0,-300.0,Interval{Float64}(-400.0,-200.0),seed2,seed2,false)
    z3 = x3*y3
    @test isapprox(z3.cc,70000,atol=1E-4)
    @test isapprox(z3.cv,50000,atol=1E-4)

    # CHECK ME AGAIN???? -47187.5 new, -47460.9375 old
    x4 = MC{2,Diff}(150.0,150.0,Interval{Float64}(100.0,200.0),seed1,seed1,false)
    y4 = MC{2,Diff}(-250.0,-250.0,Interval{Float64}(-500.0,-100.0),seed2,seed2,false)
    z4 = x4*y4
    @test isapprox(z4.cc,-30000,atol=1E-3)
    @test isapprox(z4.cv,-47187.5,atol=1E-3)

    x5 = MC{2,Diff}(-150.0,-150.0,Interval{Float64}(-200.0,-100.0),seed1,seed1,false)
    y5 = MC{2,Diff}(300.0,300.0,Interval{Float64}(200.0,400.0),seed2,seed2,false)
    z5 = x5*y5
    @test isapprox(z5.cv,-50000,atol=1E-4)
    @test isapprox(z5.cc,-40000,atol=1E-4)
end

@testset "Division" begin
    X = MC{2,NS}(3.0,3.0,Interval{Float64}(2.0,4.0), seed_gradient(1,Val(2)),seed_gradient(1,Val(2)),false)
    Y = MC{2,NS}(-4.0,-4.0,Interval{Float64}(-5.0,-3.0), seed_gradient(2,Val(2)), seed_gradient(2,Val(2)),false)
    out = X/Y
    @test isapprox(out.cc,-0.7000000000000002,atol=1E-6)
    @test isapprox(out.cv,-0.8666666666666665,atol=1E-6)
    @test isapprox(out.cc_grad[1],-0.19999999999999998,atol=1E-6)
    @test isapprox(out.cc_grad[2],-0.125,atol=1E-6)
    @test isapprox(out.cv_grad[1],-0.33333333333333337,atol=1E-6)
    @test isapprox(out.cv_grad[2],-0.1333333333333333,atol=1E-6)
    @test isapprox(out.Intv.lo,-1.33333333,atol=1E-6)
    @test isapprox(out.Intv.hi,-0.39999999999999997,atol=1E-6)
end

# REVERSE OPERATORS
#=
#=
function MC_1_is_equal(y, x, tol)
    bool1 = isapprox(y.cc,x.cc,atol=tol)
    bool2 = isapprox(y.cv,x.cv,atol=tol)
    bool3 = isapprox(y.cv_grad[1], x.cv_grad[1], atol=tol)
    bool4 = isapprox(y.cc_grad[1], x.cc_grad[1], atol=tol)
    bool5 = isapprox(y.Intv.lo, x.Intv.lo, atol=tol)
    bool6 = isapprox(y.Intv.hi, x.Intv.hi, atol=tol)
    return (bool1 && bool2 && bool3 && bool4 && bool5 && bool6)
end
=#
#=
a = MC{1}(1.0,Interval{Float64}(0.4,3.0),1)
a1 = MC{1}(-7.0,Interval{Float64}(-12.0,-4.0),1)
b = MC{1}(Interval{Float64}(-10.0,-1.0))
c = MC{1}(2.0,Interval{Float64}(1.1,4.5),1)
aout1, bout1, cout1 = mul_rev(a,b,c)
aout2, bout2, cout2 = div_rev(a,b,c)
=#

#=
@testset "Reverse Multiplication" begin

    # THE BINARY OPERATOR
    a = MC{1}(1.0, Interval{Float64}(0.4,3.0), 1)
    b = MC{1}(Interval{Float64}(-10.0,-1.0))
    c = MC{1}(2.0, Interval{Float64}(1.1,4.5), 1)

    aout1, bout1, cout1 = mul_rev(a,b,c)

    @test bout1.Intv.lo == Inf
    @test bout1.Intv.hi == -Inf
    @test cout1.Intv.lo == Inf
    @test cout1.Intv.hi == -Inf

    bout1 = MC{1}(1.0, Interval{Float64}(0.4,3.0), 1)
    cout1 = MC{1}(Interval{Float64}(-10.0,-1.0))
    aout1 = bout1*cout1

    aout1_a, bout1_a, cout1_a = mul_rev(aout1, bout1, cout1)

    MC_1_is_equal(aout1_a, aout1, 0.00001)
    MC_1_is_equal(bout1_a, bout1, 0.00001)
    MC_1_is_equal(cout1_a, cout1, 0.00001)

    bout2 = MC{1}(1.0, Interval{Float64}(0.4,3.0), 1)
    cout2 = MC{1}(Interval{Float64}(-10.0,-1.0))
    aout2 = 0.3*bout1*cout1+1.0

    aout2_a, bout2_a, cout2_a = mul_rev(aout2, bout2, cout2)

    MC_1_is_equal(aout2_a, aout2, 0.00001)
    MC_1_is_equal(bout2_a, bout2, 0.00001)
    @test cout2_a.Intv.lo == -10.0
    @test cout2_a.Intv.hi == -1.0
    @test cout2_a.cv == -6.32
    @test cout2_a.cc == -1.0
    @test cout2_a.cv_grad[1] == -8.38
    @test cout2_a.cc_grad[1] == 0.0

    # WITH FLOAT

end
=#

#=
@testset "Reverse Addition" begin
end
=#

#=
@testset "Reverse Division" begin
end
=#

#=
@testset "Reverse Subtraction" begin
end
=#

#=
@testset "Reverse Exponential" begin
    a = MC{1}(1.0,Interval{Float64}(0.4,3.0),1)
    expa = exp(a)*1.1
    y,x = exp_rev(expa,a)

    @test MC_1_is_equal(expa, y, 0.00001)
    @test x.cc == 1.0
    @test x.cv == 1.0
    @test x.Intv.lo == 0.49531
    @test x.Intv.hi == 3.0
    @test x.cc_grad == 1.0
    @test x.cv_grad == 1.0

    exp2a = exp2(a)*1.1
    y,x = exp2_rev(exp2a,a)

    @test MC_1_is_equal(exp2a, y, 0.00001)
    @test x.cc == 1.0
    @test x.cv == 1.0
    @test x.Intv.lo == 0.53753
    @test x.Intv.hi == 3.0
    @test x.cc_grad == 1.0
    @test x.cv_grad == 1.0

    exp10a = exp10(a)*1.1
    y,x = exp10_rev(exp10a,a)

    @test MC_1_is_equal(exp10a, y, 0.00001)
    @test x.cc == 1.0
    @test x.cv == 1.0
    @test x.Intv.lo == 0.441392
    @test x.Intv.hi == 3.0
    @test x.cc_grad == 1.0
    @test x.cv_grad == 1.0

    expm1a = expm1(a)*1.1
    y,x = expm1_rev(expm1a,a)

    @test MC_1_is_equal(expm1a, y, 0.00001)
    @test x.cc == 1.0
    @test x.cv == 1.0
    @test x.Intv.lo == 0.432436
    @test x.Intv.hi == 3.0
    @test x.cc_grad == 1.0
    @test x.cv_grad == 1.0
end
=#

#=
@testset "Reverse Logarithm" begin

    a = MC{1}(9.6,Interval{Float64}(9.4,10.0),1)
    a1 = MC{1}(1.0,Interval{Float64}(0.2,5.0),1)
    loga = log(a)*5.1
    y,x = log_rev(loga,a)

    @test x.Intv.lo == Inf
    @test x.Intv.hi == -Inf

end
=#

#=
# BROKEN FLOAT REVERSE
bout1 = MC{1}(1.0,Interval{Float64}(0.4,3.0),1)
cout1 = -3.0
aout1 = bout1*cout1
aout1_a, bout1_a, cout1_a = mul_rev(aout1,bout1,cout1)
=#

#=
a = MC{1}(9.6,Interval{Float64}(9.4,10.0),1)
a1 = MC{1}(1.0,Interval{Float64}(0.2,5.0),1)
loga = log(a)*5.1
y,x = log_rev(loga,a)
=#

#a0 = MC{1}(7.0,Interval{Float64}(4.5,12.0),1)
#b0 = MC{1}(Interval{Float64}(6.0,9.0))
#c0 = MC{1}(5.0,Interval{Float64}(4.1,9.5),1)
#aout3, bout3, cout3 = pow_rev(a0,b0,c0)

=#

# IMPLICIT ROUTINES
#=
@testset "Test Base Implicit Routines (Out-of-place)" begin
    EAGO.set_mc_differentiability!(0)
    opts1 = mc_opts(0.5,1,:Dense,:Newton,1,1,1E-10)

    f(x,p) = x[1]*p[1]+p[1]
    g(x,p) = [x[1]*p[1]+p[1];
              x[1]*p[1]+2*p[1]]
    function h1(x,p)
        t1 = x[1]^2
        t2 = x[1]*p[1]
        t3 = 4.0
        t4 = t1 + t2
        t5 = t4 + t3
        return [t5]
    end

    hj1(x,p) = [2*x[1]+p[1]]

    P = [Interval{Float64}(6.0,9.0)]
    X = [Interval{Float64}(-0.78,-0.4)]
    p = [7.5]
    pmid = mid.(P)

    xIntv1 = Interval{Float64}(1.0,3.0)
    xIBox = SVector{1,Interval{Float64}}([xIntv1])
    mBox = mid.(xIBox)
    np = 1
    szero = @SVector zeros(np)
    sone = @SVector ones(np)
    p_mc = [MC{np}(p[i],p[i],Interval{Float64}(P[i].lo,P[i].hi),sone,sone,false) for i=1:np]

    param = gen_expansion_params(h1,hj1,X,P,pmid,opts1)

    hbnds = implicit_relax_h(h1,hj1,p_mc,pmid,X,P,opts1,param)

    fbnds = implicit_relax_f(f,h1,hj1,X,P,p,pmid,opts1,param)

    fgbnds = implicit_relax_fg(f,g,h1,hj1,X,P,p,pmid,opts1,param)

    param = gen_expansion_params(h1,hj1,X,P,pmid,opts1)

    @test isapprox(param[2][1].cc,-0.4,atol=1E-4)
    @test isapprox(param[2][1].cv,-0.78,atol=1E-4)
    hbnds = implicit_relax_h(h1,hj1,p_mc,pmid,X,P,opts1,param)
    @test isapprox(hbnds[1].cc,-0.5209127114919797,atol=1E-4)
    @test isapprox(hbnds[1].cv,-0.6307841683146562,atol=1E-4)
    @test isapprox(hbnds[1].cc_grad[1],0.0455312,atol=1E-4)
    @test isapprox(hbnds[1].cv_grad[1],0.0941469,atol=1E-4)

    fbnds = implicit_relax_f(f,h1,hj1,X,P,p,pmid,opts1,param)
    @test isapprox(fbnds.cc,3.774523731048122,atol=1E-4)
    @test isapprox(fbnds.cv,2.5572882333553046,atol=1E-4)
    @test isapprox(fbnds.cc_grad[1],0.873187,atol=1E-4)
    @test isapprox(fbnds.cv_grad[1],0.792877,atol=1E-4)

    fgbnds = implicit_relax_fg(f,g,h1,hj1,X,P,p,pmid,opts1,param)
    @test isapprox(fgbnds[1].cc,3.774523731048122,atol=1E-4)
    @test isapprox(fgbnds[1].cv,2.5572882333553046,atol=1E-4)
    @test isapprox(fgbnds[2][1].cc,3.774523731048122,atol=1E-4)
    @test isapprox(fgbnds[2][1].cv,2.5572882333553046,atol=1E-4)
    @test isapprox(fgbnds[2][2].cc,11.274523731048122,atol=1E-4)
    @test isapprox(fgbnds[2][2].cv,10.057288233355305,atol=1E-4)

    @test isapprox(fgbnds[1].cc_grad[1],0.873187,atol=1E-4)
    @test isapprox(fgbnds[1].cv_grad[1],0.792877,atol=1E-4)
    @test isapprox(fgbnds[2][1].cc_grad[1],0.873187,atol=1E-4)
    @test isapprox(fgbnds[2][1].cv_grad[1],0.792877,atol=1E-4)
    @test isapprox(fgbnds[2][2].cc_grad[1],1.87319,atol=1E-4)
    @test isapprox(fgbnds[2][2].cv_grad[1],1.79288,atol=1E-4)
end

@testset "Test Base Implicit Routines (In-place, Single Step)" begin
    EAGO.set_mc_differentiability!(0)
    mc_opts1 = mc_opts(0.5,1,:Dense,:Newton,1,1,1E-10)

    P = [Interval{Float64}(6.0,9.0)]
    X = [Interval{Float64}(-0.78,-0.4)]
    pmid = mid.(P)
    pref_mc = [MC{1}(pmid[1],P[1],1)]
    p_mc = [MC{1}(pmid[1],P[1],1)]

    flt_param = [2.0; 1.0]

    Y = [0.0 0.0; 0.0 0.0]
    H = fill(zero(MC{1}), (1,))
    J = fill(zero(MC{1}), (2,2))

    x_mc = fill(zero(MC{1}), (1,))
    xp_mc = fill(zero(MC{1}), (1,))
    xa_mc = fill(zero(MC{1}), (1,))
    xA_mc = fill(zero(MC{1}), (1,))
    z_mc = fill(zero(MC{1}), (1,))
    aff_mc = fill(zero(MC{1}), (1,))

    interval_bnds = true
    precond = true

    z_mc =[MC{1}(Interval{Float64}(mid(X[1])))]
    t1 = z_mc[1]^2
    t2 = z_mc[1]*p_mc[1]
    t3 = 4.0
    t4 = t1 + t2 + t3

    function h3!(out, z, xp, p, t)
        out[1] = z[1]^2 + z[1]*p[1] + 4.0
    end
    function hj3!(out, z, xp, p, t)
        out[1] = 2.0*z[1]+p[1]
    end

    param = fill(zero(MC{1}), (2,))# probably this

    gen_expansion_params!(h3!, hj3!, pref_mc, xp_mc, x_mc, xa_mc, xA_mc, z_mc, aff_mc,
                          X, P, mc_opts1, param, H, J, Y, interval_bnds, flt_param, precond)
    x_mc_inplace = deepcopy(x_mc)
    params_inplace = deepcopy(param)

    implicit_relax_h!(h3!, hj3!, p_mc, pref_mc, xp_mc, x_mc, xa_mc, xA_mc, z_mc, aff_mc,
                      X, P, mc_opts1, param, H, J, Y, interval_bnds, flt_param, precond)
    hrelax_inplace = deepcopy(x_mc)

    hbnds = x_mc
    @test isapprox(hrelax_inplace[1].cc,-0.5209127114919797,atol=1E-4)
    @test isapprox(hrelax_inplace[1].cv,-0.6307841683146562,atol=1E-4)
    @test isapprox(hrelax_inplace[1].cc_grad[1],0.0455312,atol=1E-4)
    @test isapprox(hrelax_inplace[1].cv_grad[1],0.0941469,atol=1E-4)

    @test isapprox(params_inplace[2][1].cc,-0.4,atol=1E-4)
    @test isapprox(params_inplace[2][1].cv,-0.78,atol=1E-4)
end

@testset "Test Base Implicit Routines (In-place, Two Step)" begin
    EAGO.set_mc_differentiability!(0)
    mc_opts1 = mc_opts(0.5,2,:Dense,:Newton,1,1,1E-10)

    f(x, p) = x[1]*p[1]+p[1]
    g(x, p) = [x[1]*p[1]+p[1];
               x[1]*p[1]+2*p[1]]

    P = [Interval{Float64}(6.0,9.0)]
    X = [Interval{Float64}(-0.78,-0.4)]
    pmid = mid.(P)
    pref_mc = [MC{1}(pmid[1],P[1],1)]
    p_mc = [MC{1}(pmid[1],P[1],1)]

    flt_param = [2.0; 1.0]

    Y = [0.0 0.0; 0.0 0.0]
    H = fill(zero(MC{1}), (1,))
    J = fill(zero(MC{1}), (2,2))

    x_mc = fill(zero(MC{1}), (1,))
    xp_mc = fill(zero(MC{1}), (1,))
    xa_mc = fill(zero(MC{1}), (1,))
    xA_mc = fill(zero(MC{1}), (1,))
    z_mc = fill(zero(MC{1}), (1,))
    aff_mc = fill(zero(MC{1}), (1,))

    interval_bnds = true
    precond = true

    z_mc =[MC{1}(Interval{Float64}(mid(X[1])))]
    t1 = z_mc[1]^2
    t2 = z_mc[1]*p_mc[1]
    t3 = 4.0
    t4 = t1 + t2 + t3

    function h3!(out, z, xp, p, t)
        out[1] = z[1]^2 + z[1]*p[1] + 4.0
    end
    function hj3!(out, z, xp, p, t)
        out[1] = 2.0*z[1]+p[1]
    end

    param = fill(zero(MC{1}), (1,3))

    gen_expansion_params!(h3!, hj3!, pref_mc, xp_mc, x_mc, xa_mc, xA_mc, z_mc, aff_mc,
                          X, P, mc_opts1, param, H, J, Y, interval_bnds, flt_param, precond)
    x_mc_inplace = deepcopy(x_mc)
    params_inplace = deepcopy(param)

    implicit_relax_h!(h3!, hj3!, p_mc, pref_mc, xp_mc, x_mc, xa_mc, xA_mc, z_mc, aff_mc,
                      X, P, mc_opts1, param, H, J, Y, interval_bnds, flt_param, precond)
    hrelax_inplace = deepcopy(x_mc)

    hbnds = x_mc
    @test isapprox(hrelax_inplace[1].cc,-0.5209127114919797,atol=1E-4)
    @test isapprox(hrelax_inplace[1].cv,-0.6307841683146562,atol=1E-4)
    @test isapprox(hrelax_inplace[1].cc_grad[1],0.045531201933640894,atol=1E-4)
    @test isapprox(hrelax_inplace[1].cv_grad[1],0.09414689079323225,atol=1E-4)

    @test isapprox(params_inplace[2][1].cc,-0.4,atol=1E-4)
    @test isapprox(params_inplace[2][1].cv,-0.78,atol=1E-4)
end
=#
