# Create functions for comparing MC object to reference object and detail failures
function chk_ref_kernel(y::MC, yref::MC, mctol::Float64)
    pass_flag::Bool = true
    descr = "Failing Components: ("
    ~isapprox(y.cv, yref.cv; atol = mctol) && (descr = descr*" cv = $(y.cv)"; pass_flag = false)
    ~isapprox(y.cc, yref.cc; atol = mctol) && (descr = descr*" cc = $(y.cc) "; pass_flag = false)
    ~isapprox(y.Intv.lo, yref.Intv.lo; atol = mctol) && (descr = descr*" Intv.lo = $(y.Intv.lo) "; pass_flag = false)
    ~isapprox(y.Intv.hi, yref.Intv.hi; atol = mctol) && (descr = descr*" Intv.hi = $(y.Intv.hi) "; pass_flag = false)
    ~isapprox(y.cv_grad[1], yref.cv_grad[1]; atol = mctol) && (descr = descr*" cv_grad[1] "; pass_flag = false)
    ~isapprox(y.cv_grad[2], yref.cv_grad[2]; atol = mctol) && (descr = descr*" cv_grad[2] "; pass_flag = false)
    ~isapprox(y.cc_grad[1], yref.cc_grad[1]; atol = mctol) && (descr = descr*" cc_grad[1] "; pass_flag = false)
    ~isapprox(y.cc_grad[2], yref.cc_grad[2]; atol = mctol) && (descr = descr*" cc_grad[2] "; pass_flag = false)
    println(descr*")")
    pass_flag
end
check_vs_ref1(f::Function, x::MC, yref::MC, mctol::Float64) = chk_ref_kernel(f(x), yref, mctol)
check_vs_ref2(f::Function, x::MC, y::MC, yref::MC, mctol::Float64) = chk_ref_kernel(f(x,y), yref, mctol)
check_vs_refv(f::Function, x::MC, c::Float64, yref::MC, mctol::Float64) = chk_ref_kernel(f(x, c), yref, mctol)
function check_vs_ref2(f::Function, x::MC, c::Float64, yref1::MC, yref2::MC, mctol::Float64)
    pass_flag = check_vs_ref_kernel(f(x, c), yref1, mctol)
    check_vs_ref_kernel(f(c, x), yref2, mctol) && pass_flag
end

@testset "Test Univariate" begin

   mctol = 1E-4

   ##### Exponent #####
   x_exp = MC{2}(2.0, 3.0, EAGO.IntervalType(1.0,4.0), seed_gradient(Float64,1,2), seed_gradient(Float64,1,2), false)
   x_exp2 = MC{2}(2.0, 3.0, EAGO.IntervalType(1.0,4.0), seed_gradient(Float64,1,2), seed_gradient(Float64,1,2), false)
   x_exp10 = MC{2}(2.0, 3.0, EAGO.IntervalType(1.0,4.0), seed_gradient(Float64,1,2), seed_gradient(Float64,1,2), false)
   x_expm1 = MC{2}(2.0, 3.0, EAGO.IntervalType(1.0,4.0), seed_gradient(Float64,1,2), seed_gradient(Float64,1,2), false)

   yref_exp = MC{2}(15.154262241479262, 37.30486063158251, EAGO.IntervalType(2.71828, 54.5982), SVector{2,Float64}([0.0, 0.0]), SVector{2,Float64}([17.2933, 0.0]), false)
   yref_exp2 = MC{2}(4.0, 11.33333333333333, EAGO.IntervalType(1.99999, 16.0001), SVector{2,Float64}([2.77259, 0.0]), SVector{2,Float64}([4.66667, 0.0]), false)
   yref_exp10 = MC{2}(100.0, 6670.0000000000055, EAGO.IntervalType(9.999999999999999999, 10000.00000000001), SVector{2,Float64}([2302.5850929940457, 0.0]), SVector{2,Float64}([3330.0, 0.0]), false)
   yref_expm1 = MC{2}(6.38905609893065, 36.304860631582514, EAGO.IntervalType(1.71828, 53.5982), SVector{2,Float64}([7.38906, 0.0]), SVector{2,Float64}([17.2933, 0.0]), false)

   EAGO.set_mc_differentiability!(0); @test check_vs_ref1(exp, x_exp, yref_exp, mctol)
   EAGO.set_mc_differentiability!(0); @test check_vs_ref1(exp2, x_exp2, yref_exp2, mctol)
   EAGO.set_mc_differentiability!(0); @test check_vs_ref1(exp10, x_exp10, yref_exp10, mctol)
   EAGO.set_mc_differentiability!(0); @test check_vs_ref1(expm1, x_expm1, yref_expm1, mctol)

   ##### Logarithm #####
   x_log = MC{2}(2.0, 3.0, EAGO.IntervalType(1.0,4.0), seed_gradient(Float64,1,2), seed_gradient(Float64,1,2), false)
   x_log2 = MC{2}(2.0, 3.0, EAGO.IntervalType(1.0,4.0), seed_gradient(Float64,1,2), seed_gradient(Float64,1,2), false)
   x_log10 = MC{2}(2.0, 3.0, EAGO.IntervalType(1.0,4.0), seed_gradient(Float64,1,2), seed_gradient(Float64,1,2), false)
   x_log1p = MC{2}(2.0, 3.0, EAGO.IntervalType(1.0,4.0), seed_gradient(Float64,1,2), seed_gradient(Float64,1,2), false)

   yref_log = MC{2}(0.46209812037329695, 0.6931471805599453, EAGO.IntervalType(0, 1.3863), SVector{2,Float64}([0.462098, 0.0]), SVector{2,Float64}([0.5, 0.0]), false)
   yref_log2 = MC{2}(0.6666666666666666, 1.0, EAGO.IntervalType(0, 2), SVector{2,Float64}([0.666667, 0.0]), SVector{2,Float64}([0.721348, 0.0]), false)
   yref_log10 = MC{2}(0.20068666377598746, 0.3010299956639812, EAGO.IntervalType(0, 0.60206), SVector{2,Float64}([0.200687, 0.0]), SVector{2,Float64}([0.217147, 0.0]), false)
   yref_log1p = MC{2}(0.998577424517997, 1.0986122886681098, EAGO.IntervalType(0.693147, 1.60944), SVector{2,Float64}([0.30543, 0.0]), SVector{2,Float64}([0.333333, 0.0]), false)

   EAGO.set_mc_differentiability!(0); @test check_vs_ref1(log, x_log, yref_log, mctol)
   EAGO.set_mc_differentiability!(0); @test check_vs_ref1(log2, x_log2, yref_log2, mctol)
   EAGO.set_mc_differentiability!(0); @test check_vs_ref1(log10, x_log10, yref_log10, mctol)
   EAGO.set_mc_differentiability!(0); @test check_vs_ref1(log1p, x_log1p, yref_log1p, mctol)

   #####  Square root #####
   x_sqrt_ns = MC{2}(4.5, 4.5, EAGO.IntervalType(3.0,9.0), seed_gradient(Float64,1,2), seed_gradient(Float64,1,2), false)
   x_sqrt_d1 = MC{2}(4.0, 4.0, EAGO.IntervalType(3.0,7.0), seed_gradient(Float64,1,2), seed_gradient(Float64,1,2), false)

   yref_sqrt_ns = MC{2}(2.049038105676658, 2.1213203435596424, EAGO.IntervalType(1.73205, 3.0), SVector{2,Float64}([0.211325, 0.0]), SVector{2,Float64}([0.235702, 0.0]), false)
   yref_sqrt_d1 = MC{2}(1.9604759334428057, 2.0, EAGO.IntervalType(1.73205, 2.64576), SVector{2,Float64}([0.228425, 0.0]), SVector{2,Float64}([0.25, 0.0]), false)

   EAGO.set_mc_differentiability!(0); @test check_vs_ref1(sqrt, x_sqrt_ns, yref_sqrt_ns, mctol)
   EAGO.set_mc_differentiability!(1); @test check_vs_ref1(sqrt, x_sqrt_d1, yref_sqrt_d1, mctol)

   #####  Step #####
   x_step_p = MC{2}(4.0, 4.0, EAGO.IntervalType(3.0,7.0), seed_gradient(Float64,1,2), seed_gradient(Float64,1,2), false)
   x_step_n = MC{2}(-4.0, -4.0, -EAGO.IntervalType(3.0,7.0), seed_gradient(Float64,1,2), seed_gradient(Float64,1,2), false)
   x_step_z = MC{2}(-2.0, -2.0, EAGO.IntervalType(-3.0,1.0), seed_gradient(Float64,1,2), seed_gradient(Float64,1,2), false)
   x_step_zp = MC{2}(0.5, 0.5, EAGO.IntervalType(-3.0,1.0), seed_gradient(Float64,1,2), seed_gradient(Float64,1,2), false)

   yref_step_p = MC{2}(1.0, 1.0, EAGO.IntervalType(1.0, 1.0), @SVector[0.0, 0.0], @SVector[0.0, 0.0], false)
   yref_step_n = MC{2}(0.0, 0.0, EAGO.IntervalType(0.0, 0.0), @SVector[0.0, 0.0], @SVector[0.0, 0.0], false)
   yref_step_ns_z = MC{2}(0.0, 0.3333333333333333, EAGO.IntervalType(0.0, 1.0), @SVector[0.0, 0.0], @SVector[-0.6666666666666666, 0.0], false)
   yref_step_ns_zp = MC{2}(0.5, 1.0, EAGO.IntervalType(0.0, 1.0), @SVector[1.0, 0.0], @SVector[0.0, 0.0], false)
   yref_step_d1_z = MC{2}(0.0, 0.5555555555555556, EAGO.IntervalType(0.0, 1.0), @SVector[0.0, 0.0], @SVector[0.444444, 0.0], false)

   EAGO.set_mc_differentiability!(0); @test check_vs_ref1(step, x_step_p, yref_step_p, mctol)
   EAGO.set_mc_differentiability!(0); @test check_vs_ref1(step, x_step_n, yref_step_n, mctol)
   EAGO.set_mc_differentiability!(0); @test check_vs_ref1(step, x_step_z, yref_step_ns_z, mctol)
   EAGO.set_mc_differentiability!(0); @test check_vs_ref1(step, x_step_zp, yref_step_ns_zp, mctol)
   EAGO.set_mc_differentiability!(1); @test check_vs_ref1(step, x_step_p, yref_step_p, mctol)
   EAGO.set_mc_differentiability!(1); @test check_vs_ref1(step, x_step_n, yref_step_n, mctol)
   EAGO.set_mc_differentiability!(1); @test check_vs_ref1(step, x_step_z, yref_step_d1_z, mctol)

   #####  Absolute value #####
   x_abs_ns = MC{2}(4.5, 4.5, EAGO.IntervalType(-3.0,8.0), seed_gradient(Float64,1,2), seed_gradient(Float64,1,2), false)
   x_abs_d1 = MC{2}(4.0, 4.0, EAGO.IntervalType(3.0,7.0), seed_gradient(Float64,1,2), seed_gradient(Float64,1,2), false)

   yref_abs_ns = MC{2}(4.5, 6.409090909090908, EAGO.IntervalType(0.0, 8.0), @SVector[1.0, 0.0], @SVector[0.454545, 0.0], false)
   yref_abs_d1 = MC{2}(1.3061224489795915, 2.0, EAGO.IntervalType(3.0, 7.0), @SVector[0.979592, 0.0], @SVector[1.0, 0.0], false)

   EAGO.set_mc_differentiability!(0); @test check_vs_ref1(abs, x_abs_ns, yref_abs_ns, mctol)
   EAGO.set_mc_differentiability!(1); @test check_vs_ref1(abs, x_abs_d1, yref_abs_d1, mctol)

   #####  Sign #####
   x_sign_d1_p = MC{2}(4.0, 4.0, EAGO.IntervalType(3.0,7.0), seed_gradient(Float64,1,2), seed_gradient(Float64,1,2), false)
   x_sign_d1_n = MC{2}(-4.0, -4.0, EAGO.IntervalType(-7.0,-3.0), seed_gradient(Float64,1,2), seed_gradient(Float64,1,2), false)
   x_sign_d1_z = MC{2}(-2.0, -2.0, EAGO.IntervalType(-3.0,1.0), seed_gradient(Float64,1,2), seed_gradient(Float64,1,2), false)

   yref_sign_d1_p = MC{2}(1.0, 1.0, EAGO.IntervalType(1.0, 1.0), @SVector[0.0, 0.0], @SVector[0.0, 0.0], false)
   yref_sign_d1_n = MC{2}(-1.0, -1.0, EAGO.IntervalType(-1.0, -1.0), @SVector[0.0, 0.0], @SVector[0.0, 0.0], false)
   yref_sign_d1_z = MC{2}(-1.0, 0.11111111111111116, EAGO.IntervalType(-1.0, 1.0), @SVector[0.0, 0.0], @SVector[0.888889, 0.0], false)

   EAGO.set_mc_differentiability!(1); @test check_vs_ref1(sign, x_sign_d1_p, yref_sign_d1_p, mctol)
   EAGO.set_mc_differentiability!(1); @test check_vs_ref1(sign, x_sign_d1_n, yref_sign_d1_n, mctol)
   EAGO.set_mc_differentiability!(1); @test check_vs_ref1(sign, x_sign_d1_z, yref_sign_d1_z, mctol)

   #####  Sine #####
   x_sin_p = MC{2}(4.0, 4.0, EAGO.IntervalType(3.0,7.0), seed_gradient(Float64,1,2), seed_gradient(Float64,1,2), false)
   x_sin_n = MC{2}(-4.0, -4.0, EAGO.IntervalType(-7.0,-3.0), seed_gradient(Float64,1,2), seed_gradient(Float64,1,2), false)
   x_sin_z = MC{2}(-2.0, -2.0, EAGO.IntervalType(-3.0,1.0), seed_gradient(Float64,1,2), seed_gradient(Float64,1,2), false)

   yref_sin_ns_p = MC{2}(-0.7568024953079283, 0.2700866557245978, EAGO.IntervalType(-1.0, 0.656987), @SVector[-0.653644, 0.0], @SVector[0.128967, 0.0], false)
   yref_sin_ns_n = MC{2}(-0.2700866557245979, 0.7568024953079283, EAGO.IntervalType(-0.656987, 1.0), @SVector[0.128967, 0.0], @SVector[-0.653644, 0.0], false)
   yref_sin_ns_z = MC{2}(-0.9092974268256817, 0.10452774015707458, EAGO.IntervalType(-1.0, 0.841471), @SVector[-0.416147, 0.0], @SVector[0.245648, 0.0], false)
   yref_sin_d1_z = MC{2}(-0.9092974268256817, 0.10452774015707458, EAGO.IntervalType(-1.0, 0.841471), @SVector[-0.416147, 0.0], @SVector[0.245648, 0.0], false)

   EAGO.set_mc_differentiability!(1); @test check_vs_ref1(sin, x_sin_p, yref_sin_ns_p, mctol)
   EAGO.set_mc_differentiability!(1); @test check_vs_ref1(sin, x_sin_n, yref_sin_ns_n, mctol)
   EAGO.set_mc_differentiability!(1); @test check_vs_ref1(sin, x_sin_z, yref_sin_ns_z, mctol)
   EAGO.set_mc_differentiability!(0); @test check_vs_ref1(sin, x_sin_z, yref_sin_d1_z, mctol)

   #####  Cosine #####
   x_cos_p = MC{2}(4.0, 4.0, EAGO.IntervalType(3.0,7.0), seed_gradient(Float64,1,2), seed_gradient(Float64,1,2), false)
   x_cos_n = MC{2}(-4.0, -4.0, EAGO.IntervalType(-7.0,-3.0), seed_gradient(Float64,1,2), seed_gradient(Float64,1,2), false)
   x_cos_z = MC{2}(-2.0, -2.0, EAGO.IntervalType(-3.0,1.0), seed_gradient(Float64,1,2), seed_gradient(Float64,1,2), false)

   yref_cos_d1_p = MC{2}(-0.703492113936536, -0.31034065427934965, EAGO.IntervalType(-1.0, 1.0), @SVector[0.485798, 0.0], @SVector[0.679652, 0.0], false)
   yref_cos_d1_n = MC{2}(-0.703492113936536, -0.31034065427934965, EAGO.IntervalType(-1.0, 1.0), @SVector[-0.485798, 0.0], @SVector[-0.679652, 0.0], false)
   yref_cos_d1_z = MC{2}(-0.6314158569813042, -0.222468094224762, EAGO.IntervalType(-0.989993, 1.0), @SVector[0.390573, 0.0], @SVector[0.76752, 0.0], false)
   yref_cos_ns_z = MC{2}(-0.6314158569813042, -0.222468094224762, EAGO.IntervalType(-0.989993, 1.0), @SVector[0.390573, 0.0], @SVector[0.76752, 0.0], false)

   EAGO.set_mc_differentiability!(1); @test check_vs_ref1(cos, x_cos_p, yref_cos_d1_p, mctol)
   EAGO.set_mc_differentiability!(1); @test check_vs_ref1(cos, x_cos_n, yref_cos_d1_n, mctol)
   EAGO.set_mc_differentiability!(1); @test check_vs_ref1(cos, x_cos_z, yref_cos_d1_z, mctol)
   EAGO.set_mc_differentiability!(0); @test check_vs_ref1(cos, x_cos_z, yref_cos_ns_z, mctol)

   ##### Tangent #####
   x_tan_p = MC{2}(0.6, 0.6, EAGO.IntervalType(0.5,1.0), seed_gradient(Float64,1,2), seed_gradient(Float64,1,2), false)
   x_tan_n = MC{2}(-0.8, -0.8, EAGO.IntervalType(-1.0,-0.5), seed_gradient(Float64,1,2), seed_gradient(Float64,1,2), false)
   x_tan_z = MC{2}(-0.3,-0.3, EAGO.IntervalType(-0.5,0.5), seed_gradient(Float64,1,2), seed_gradient(Float64,1,2), false)
   x_tan_err = MC{2}(0.6, 0.6, EAGO.IntervalType(-4.5, 5.0), seed_gradient(Float64,1,2), seed_gradient(Float64,1,2), false)

   yref_tan_ns_p = MC{2}(0.6841368083416923, 0.7485235368060128, EAGO.IntervalType(0.546302, 1.55741), @SVector[1.46804, 0.0], @SVector[2.02221, 0.0], false)
   yref_tan_ns_n = MC{2}(-1.1529656307304577, -1.0296385570503641, EAGO.IntervalType(-1.55741, -0.546302), @SVector[2.02221, 0.0], @SVector[2.06016, 0.0], false)
   yref_tan_ns_z = MC{2}(-0.332534, -0.30933624960962325, EAGO.IntervalType(-0.546303,0.546303), @SVector[1.06884, 0.0], @SVector[1.09569, 0.0], false)
   yref_tan_d1_z = MC{2}(-0.332534, -0.309336, EAGO.IntervalType(-0.546303,0.546303), @SVector[1.06884, 0.0], @SVector[1.09569, 0.0], false)

   EAGO.set_mc_differentiability!(0); @test check_vs_ref1(tan, x_tan_p, yref_tan_ns_p, mctol)
   EAGO.set_mc_differentiability!(0); @test check_vs_ref1(tan, x_tan_n, yref_tan_ns_n, mctol)
   EAGO.set_mc_differentiability!(0); @test check_vs_ref1(tan, x_tan_z, yref_tan_ns_z, mctol)
   EAGO.set_mc_differentiability!(1); @test check_vs_ref1(tan, x_tan_z, yref_tan_d1_z, mctol)
   @test_throws ErrorException tan(x_tan_err)

   #####  Arcsine #####
   x_asin_p = MC{2}(-0.7, -0.7, EAGO.IntervalType(-0.9,-0.5), seed_gradient(Float64,1,2), seed_gradient(Float64,1,2), false)
   x_asin_n = MC{2}(0.7, 0.7, EAGO.IntervalType(0.5,0.9), seed_gradient(Float64,1,2), seed_gradient(Float64,1,2), false)
   x_asin_z = MC{2}(-0.1, -0.1, EAGO.IntervalType(-0.5,0.5), seed_gradient(Float64,1,2), seed_gradient(Float64,1,2), false)

   yref_asin_d1_p = MC{2}(-0.8216841452984665, -0.775397496610753, EAGO.IntervalType(-1.11977, -0.523598), @SVector[1.49043, 0.0], @SVector[1.40028, 0.0], false)
   yref_asin_d1_n = MC{2}(0.775397496610753, 0.8216841452984665, EAGO.IntervalType(0.523598, 1.11977), @SVector[1.40028, 0.0], @SVector[1.49043, 0.0], false)
   yref_asin_d1_z = MC{2}(-0.10958805193420748, -0.0974173098978382, EAGO.IntervalType(-0.523599, 0.523599), @SVector[1.03503, 0.0], @SVector[1.03503, 0.0], false)
   yref_asin_ns_z = MC{2}(-0.10958805193420748, -0.0974173098978382, EAGO.IntervalType(-0.523599, 0.523599), @SVector[1.03503, 0.0], @SVector[1.03503, 0.0], false)

   EAGO.set_mc_differentiability!(1); @test check_vs_ref1(asin, x_asin_p, yref_asin_d1_p, mctol)
   EAGO.set_mc_differentiability!(1); @test check_vs_ref1(asin, x_asin_n, yref_asin_d1_n, mctol)
   EAGO.set_mc_differentiability!(1); @test check_vs_ref1(asin, x_asin_z, yref_asin_d1_z, mctol)
   EAGO.set_mc_differentiability!(0); @test check_vs_ref1(asin, x_asin_z, yref_asin_ns_z, mctol)

   ##### Arctangent #####
   x_atan_p = MC{2}(4.0, 4.0, EAGO.IntervalType(3.0, 7.0), seed_gradient(Float64,1,2), seed_gradient(Float64,1,2), false)
   yref_atan_d1_p = MC{2}(1.294009147346374, 1.3258176636680326, EAGO.IntervalType(1.24904, 1.4289), seed_gradient(Float64,1,2), seed_gradient(Float64,1,2), false)
   EAGO.set_mc_differentiability!(1); @test check_vs_ref1(atan, x_atan_p, yref_atan_d1_p, mctol)

   ##### Hyperbolic Sine #####
   x_sinh_p = MC{2}(4.0, 4.0, EAGO.IntervalType(3.0,7.0), seed_gradient(Float64,1,2), seed_gradient(Float64,1,2), false)
   x_sinh_n = MC{2}(4.0, -4.0, EAGO.IntervalType(-7.0,-3.0), seed_gradient(Float64,1,2), seed_gradient(Float64,1,2), false)
   x_sinh_z = MC{2}(-2.0, -2.0, EAGO.IntervalType(-3.0,1.0), seed_gradient(Float64,1,2), seed_gradient(Float64,1,2), false)

   yref_sinh_d1_p = MC{2}(27.28991719712775, 144.59243701386904, EAGO.IntervalType(10.0178, 548.3161232732466), @SVector[27.3082, 0.0], @SVector[134.575, 0.0], false)
   yref_sinh_d1_n = MC{2}(-144.59243701386904, -10.017874927409903, EAGO.IntervalType(-548.3161232732466, -10.0178), @SVector[134.575, 0.0], @SVector[27.3082, 0.0], false)
   yref_sinh_d1_z = MC{2}(-7.219605897146477, -3.626860407847019, EAGO.IntervalType(-10.0179, 1.17521), @SVector[2.79827, 0.0], @SVector[3.7622, 0.0], false)
   yref_sinh_ns_p = MC{2}(27.28991719712775, 144.59243701386904, EAGO.IntervalType(10.0178, 548.317), @SVector[27.3082, 0.0], @SVector[134.575, 0.0], false)

   EAGO.set_mc_differentiability!(1); @test check_vs_ref1(sinh, x_sinh_p, yref_sinh_d1_p, mctol)
   EAGO.set_mc_differentiability!(1); @test check_vs_ref1(sinh, x_sinh_n, yref_sinh_d1_n, mctol)
   EAGO.set_mc_differentiability!(1); @test check_vs_ref1(sinh, x_sinh_z, yref_sinh_d1_z, mctol)
   EAGO.set_mc_differentiability!(0); @test check_vs_ref1(sinh, x_sinh_p, yref_sinh_d1_p, mctol)

   ##### Hyperbolic Cosine #####
   x_cosh = MC{2}(4.0, 4.0, EAGO.IntervalType(3.0,7.0), seed_gradient(Float64,1,2), seed_gradient(Float64,1,2), false)
   yref_cosh_ns = MC{2}(27.308232836016487, 144.63000528563632, EAGO.IntervalType(10.0676, 548.318), @SVector[-27.2899, 0.0], @SVector[134.562, 0.0], false)
   yref_cosh_d1 = MC{2}(27.308232836016487, 144.63000528563632, EAGO.IntervalType(10.0676, 548.318), @SVector[-27.2899, 0.0], @SVector[134.562, 0.0], false)
   EAGO.set_mc_differentiability!(0); @test check_vs_ref1(cosh, x_cosh, yref_cosh_ns, mctol)
   EAGO.set_mc_differentiability!(1); @test check_vs_ref1(cosh, x_cosh, yref_cosh_d1, mctol)

   ##### Hyperbolic Tangent #####
   x_tanh_p = MC{2}(4.0, 4.0, EAGO.IntervalType(3.0,7.0), seed_gradient(Float64,1,2), seed_gradient(Float64,1,2), false)
   x_tanh_n = MC{2}(-4.0, -4.0, EAGO.IntervalType(-7.0,-3.0), seed_gradient(Float64,1,2), seed_gradient(Float64,1,2), false)
   x_tanh_z1 = MC{2}(-2.0, -2.0, EAGO.IntervalType(-3.0,1.0), seed_gradient(Float64,1,2), seed_gradient(Float64,1,2), false)
   x_tanh_z2 = MC{2}(2.0, 2.0, EAGO.IntervalType(-1.0,3.0), seed_gradient(Float64,1,2), seed_gradient(Float64,1,2), false)

   yref_tanh_d1_p = MC{2}(0.996290649501034, 0.999329299739067, EAGO.IntervalType(0.995054,0.999999), @SVector[0.0012359, 0.0], @SVector[0.00134095, 0.0], false)
   yref_tanh_d1_n = MC{2}(-0.999329299739067, -0.996290649501034, EAGO.IntervalType(-0.999999, -0.995054), @SVector[0.00134095, 0.0], @SVector[0.0012359, 0.0], false)
   yref_tanh_d1_z1 = MC{2}(-0.9640275800758169, -0.5558207301372651, EAGO.IntervalType(-0.995055, 0.761595), @SVector[0.0706508, 0.0], @SVector[0.439234, 0.0], false)
   yref_tanh_d1_z2 = MC{2}(0.5558207301372651, 0.9640275800758169, EAGO.IntervalType(-0.761595, 0.995055), @SVector[0.439234, 0.0], @SVector[0.0706508, 0.0], false)
   yref_tanh_ns = MC{2}(-0.9640275800758169, -0.5558207301372651, EAGO.IntervalType(-0.995055, 0.761595), @SVector[0.0706508, 0.0], @SVector[0.439234, 0.0], false)

   EAGO.set_mc_differentiability!(1); @test check_vs_ref1(tanh, x_tanh_p, yref_tanh_d1_p, mctol)
   EAGO.set_mc_differentiability!(1); @test check_vs_ref1(tanh, x_tanh_n, yref_tanh_d1_n, mctol)
   EAGO.set_mc_differentiability!(1); @test check_vs_ref1(tanh, x_tanh_z1, yref_tanh_d1_z1, mctol)
   EAGO.set_mc_differentiability!(1); @test check_vs_ref1(tanh, x_tanh_z2, yref_tanh_d1_z2, mctol)
   EAGO.set_mc_differentiability!(0); @test check_vs_ref1(tanh, x_tanh_z1, yref_tanh_ns, mctol)

   ##### Inverse Hyperbolic Sine #####
   x_asinh_p = MC{2}(0.3, 0.3, EAGO.IntervalType(0.1,0.7), seed_gradient(Float64,1,2), seed_gradient(Float64,1,2), false)
   x_asinh_n = MC{2}(-0.3, -0.3, EAGO.IntervalType(-0.7,-0.1), seed_gradient(Float64,1,2), seed_gradient(Float64,1,2), false)
   x_asinh_z1 = MC{2}(2.0, 2.0, EAGO.IntervalType(-3.0,7.0), seed_gradient(Float64,1,2), seed_gradient(Float64,1,2), false)
   x_asinh_z2 = MC{2}(-2.0, -2.0, EAGO.IntervalType(-7.0,3.0), seed_gradient(Float64,1,2), seed_gradient(Float64,1,2), false)

   yref_asinh_d1_p = MC{2}(0.2841115746269236, 0.29567304756342244, EAGO.IntervalType(0.099834,0.652667), @SVector[0.921387, 0.0], @SVector[0.957826, 0.0], false)
   yref_asinh_d1_n = MC{2}(-0.29567304756342244, -0.2841115746269236, EAGO.IntervalType(-0.652667,-0.099834), @SVector[0.957826, 0.0], @SVector[0.921387, 0.0], false)
   yref_asinh_d1_z1 = MC{2}(0.3730697449603356, 1.4436354751788103, EAGO.IntervalType(-1.81845, 2.64413), @SVector[0.45421, 0.0], @SVector[0.447214, 0.0], false)
   yref_asinh_d1_z2 = MC{2}(-1.4436354751788103, -0.3730697449603356, EAGO.IntervalType(-2.64413,1.81845), @SVector[0.447214, 0.0], @SVector[0.45421, 0.0], false)
   yref_asinh_ns = MC{2}(0.2841115746269236, 0.29567304756342244, EAGO.IntervalType(0.099834,0.652667), @SVector[0.921387, 0.0], @SVector[0.957826, 0.0], false)

   EAGO.set_mc_differentiability!(1); @test check_vs_ref1(asinh, x_asinh_p, yref_asinh_d1_p, mctol)
   EAGO.set_mc_differentiability!(1); @test check_vs_ref1(asinh, x_asinh_n, yref_asinh_d1_n, mctol)
   EAGO.set_mc_differentiability!(1); @test check_vs_ref1(asinh, x_asinh_z1, yref_asinh_d1_z1, mctol)
   EAGO.set_mc_differentiability!(1); @test check_vs_ref1(asinh, x_asinh_z2, yref_asinh_d1_z2, mctol)
   EAGO.set_mc_differentiability!(0); @test check_vs_ref1(asinh, x_asinh_z1, yref_asinh_ns, mctol)

   ##### Inverse Hyperbolic Cosine #####
   x_acosh = MC{2}(4.0, 4.0, EAGO.IntervalType(3.0,7.0), seed_gradient(Float64,1,2), seed_gradient(Float64,1,2), false)
   yref_acosh = MC{2}(1.9805393289917226, 2.0634370688955608, EAGO.IntervalType(1.76274,2.63392), @SVector[0.217792, 0.0], @SVector[0.258199, 0.0], false)
   EAGO.set_mc_differentiability!(0); @test check_vs_ref1(acosh, x_acosh, yref_acosh, mctol)

   ##### Inverse Hyperbolic Tangent #####
   x_atanh_p = MC{2}(0.3, 0.3, EAGO.IntervalType(0.1,0.7), seed_gradient(Float64,1,2), seed_gradient(Float64,1,2), false)
   x_atanh_n = MC{2}(-0.3, -0.3, EAGO.IntervalType(-0.7,-0.1), seed_gradient(Float64,1,2), seed_gradient(Float64,1,2), false)
   x_atanh_z1 = MC{2}(0.2, 0.2, EAGO.IntervalType(-0.6,0.7), seed_gradient(Float64,1,2), seed_gradient(Float64,1,2), false)
   x_atanh_z2 = MC{2}(-0.2, -0.2, EAGO.IntervalType(-0.7,0.6), seed_gradient(Float64,1,2), seed_gradient(Float64,1,2), false)

   yref_atanh_d1_p = MC{2}(0.30951960420311175, 0.3559904077187347, EAGO.IntervalType(0.100335,0.867301), @SVector[1.0989, 0.0], @SVector[1.27828, 0.0], false)
   yref_atanh_d1_n = MC{2}(-0.3559904077187347, -0.30951960420311175, EAGO.IntervalType(-0.867301,-0.100335), @SVector[1.27828, 0.0], @SVector[1.0989, 0.0], false)
   yref_atanh_d1_z1 = MC{2}(0.2027325540540822, 0.2788904617454707, EAGO.IntervalType(-0.30952,0.867301), @SVector[1.04167, 0.0], @SVector[1.17682, 0.0], false)
   yref_atanh_d1_z2 = MC{2}(-0.2788904617454707, -0.2027325540540822, EAGO.IntervalType(-0.867301,0.30952), @SVector[1.17682, 0.0], @SVector[1.04167, 0.0], false)
   yref_atanh_ns = MC{2}(0.30951960420311175, 0.3559904077187347, EAGO.IntervalType(0.100335,0.867301), @SVector[1.0989, 0.0], @SVector[1.27828, 0.0], false)

   EAGO.set_mc_differentiability!(1); @test check_vs_ref1(atanh, x_atanh_p, yref_atanh_d1_p, mctol)
   EAGO.set_mc_differentiability!(1); @test check_vs_ref1(atanh, x_atanh_n, yref_atanh_d1_n, mctol)
   EAGO.set_mc_differentiability!(1); @test check_vs_ref1(atanh, x_atanh_z1, yref_atanh_d1_z1, mctol)
   EAGO.set_mc_differentiability!(1); @test check_vs_ref1(atanh, x_atanh_z2, yref_atanh_d1_z2, mctol)
   EAGO.set_mc_differentiability!(0); @test check_vs_ref1(atanh, x_atanh_z1, yref_atanh_ns, mctol)

   # CONVERSION
   X = MC{2}(4.5, 4.5, EAGO.IntervalType(-3.0,8.0), seed_gradient(Float64, 1, 2), seed_gradient(Float64, 1, 2), false)
   X1 = convert(MC{2}, 1)
   X2 = convert(MC{2}, 1.1)
   X3 = convert(MC{2}, EAGO.IntervalType(2.1,4.3))
   @test X1.cc == 1.0
   @test X1.cv == 1.0
   @test X2.cc == 1.1
   @test X2.cv == 1.1
   @test X3.cc == 4.3
   @test X3.cv == 2.1

   @test +x_step_p == x_step_p
end

@testset "Test Multivariant w/Constant" begin
   x = MC{2}(4.5, 4.5, EAGO.IntervalType(-3.0,8.0), seed_gradient(Float64,1,2), seed_gradient(Float64,1,2), false)

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

   x1 = 1.0/x
   @test isapprox(x1.cc, -0.25, atol=mctol)
   @test isapprox(x1.cv, -0.266666666, atol=mctol)
   @test isapprox(x1.cc_grad[1], 0.0, atol=mctol)
   @test isapprox(x1.cc_grad[2], -0.0625,atol=mctol)
   @test isapprox(x1.cv_grad[1], 0.0, atol=mctol)
   @test isapprox(x1.cv_grad[2],-0.0666667, atol=mctol)

   x2 = 1/x
   @test isapprox(x2.cc, -0.25, atol=1E-6)
   @test isapprox(x2.cv, -0.266666666, atol=1E-6)
   @test isapprox(x2.cc_grad[1], 0.0, atol=1E-6)
   @test isapprox(x2.cc_grad[2], -0.0625, atol=1E-6)
   @test isapprox(x2.cv_grad[1], 0.0, atol=1E-6)
   @test isapprox(x2.cv_grad[2], -0.0666667, atol=1E-6)

   yref1_min_ns = MC{2}(1.0909090909090908, 3.0, EAGO.IntervalType(-3.0, 3.0), @SVector[0.545455, 0.0], @SVector[0.0, 0.0], false)
   yref2_min_ns = MC{2}(1.0909090909090908, 3.0, EAGO.IntervalType(-3.0, 3.0), @SVector[0.545455, 0.0], @SVector[0.0, 0.0], false)
   yref1_max_ns = MC{2}(5.0, 7.045454545454545, EAGO.IntervalType(5.0, 8.0), @SVector[0.0, 0.0], @SVector[0.272727, 0.0], false)
   yref2_max_ns = MC{2}(5.0, 7.045454545454545, EAGO.IntervalType(5.0, 8.0), @SVector[0.0, 0.0], @SVector[0.272727, 0.0], false)
   yref1_min_d1 = MC{2}(1.0909090909090908, 3.0, EAGO.IntervalType(-3.0, 3.0), @SVector[0.545455, 0.0], @SVector[0.0, 0.0], false)
   yref2_min_d1 = MC{2}(1.0909090909090908, 3.0, EAGO.IntervalType(-3.0, 3.0), @SVector[0.545455, 0.0], @SVector[0.0, 0.0], false)
   yref1_max_d1 = MC{2}(5.0, 7.045454545454545, EAGO.IntervalType(5.0, 8.0), @SVector[0.0, 0.0], @SVector[0.272727, 0.0], false)
   yref2_max_d1 = MC{2}(5.0, 7.045454545454545, EAGO.IntervalType(5.0, 8.0), @SVector[0.0, 0.0], @SVector[0.272727, 0.0], false)

   EAGO.set_mc_differentiability!(0); @test check_vs_ref2(min, x_min, 3.0, yref1_min_ns, yref2_min_ns, mctol)
   EAGO.set_mc_differentiability!(0); @test check_vs_ref2(max, x_max, 5.0, yref1_max_ns, yref2_max_ns, mctol)
   EAGO.set_mc_differentiability!(1); @test check_vs_ref2(min, x_min, 3.0, yref1_min_d1, yref2_min_d1, mctol)
   EAGO.set_mc_differentiability!(1); @test check_vs_ref2(max, x_max, 5.0, yref1_max_d1, yref2_max_d1, mctol)
   EAGO.set_mc_differentiability!(0); @test check_vs_ref2(min, x_min, Float16(3.0), yref1_min_ns, yref2_min_ns, mctol)
   EAGO.set_mc_differentiability!(0); @test check_vs_ref2(max, x_max, Float16(5.0), yref1_max_ns, yref2_max_ns, mctol)
   EAGO.set_mc_differentiability!(1); @test check_vs_ref2(min, x_min, Float16(3.0), yref1_min_d1, yref2_min_d1, mctol)
   EAGO.set_mc_differentiability!(1); @test check_vs_ref2(max, x_max, Float16(5.0), yref1_max_d1, yref2_max_d1, mctol)
   EAGO.set_mc_differentiability!(0); @test check_vs_ref2(min, x_min, Float32(3.0), yref1_min_ns, yref2_min_ns, mctol)
   EAGO.set_mc_differentiability!(0); @test check_vs_ref2(max, x_max, Float32(5.0), yref1_max_ns, yref2_max_ns, mctol)
   EAGO.set_mc_differentiability!(1); @test check_vs_ref2(min, x_min, Float32(3.0), yref1_min_d1, yref2_min_d1, mctol)
   EAGO.set_mc_differentiability!(1); @test check_vs_ref2(max, x_max, Float32(5.0), yref1_max_d1, yref2_max_d1, mctol)
   EAGO.set_mc_differentiability!(0); @test check_vs_ref2(min, x_min, 3, yref1_min_ns, yref2_min_ns, mctol)
   EAGO.set_mc_differentiability!(0); @test check_vs_ref2(max, x_max, 5, yref1_max_ns, yref2_max_ns, mctol)
   EAGO.set_mc_differentiability!(1); @test check_vs_ref2(min, x_min, 3, yref1_min_d1, yref2_min_d1, mctol)
   EAGO.set_mc_differentiability!(1); @test check_vs_ref2(max, x_max, 5, yref1_max_d1, yref2_max_d1, mctol)

   x1 = MC{2}(4.0, 4.0, EAGO.IntervalType(3.0, 7.0), seed_gradient(Float64,1,2), seed_gradient(Float64,1,2), false)
   x2 = MC{2}(-4.5, -4.5, EAGO.IntervalType(-8.0, -3.0), seed_gradient(Float64,1,2), seed_gradient(Float64,1,2), false)
   x3 = MC{2}(4.5, 4.5, EAGO.IntervalType(3.0, 8.0), seed_gradient(Float64,1,2), seed_gradient(Float64,1,2), false)
   x4 = MC{2}(-4.5, -4.5, EAGO.IntervalType(-8.0, -3.0), seed_gradient(Float64,1,2), seed_gradient(Float64,1,2), false)
   x5 = MC{2}(-4.0, -4.0, EAGO.IntervalType(-7.0, -3.0), seed_gradient(Float64,1,2), seed_gradient(Float64,1,2), false)
   x6 = MC{2}(-2.0, -2.0,EAGO.IntervalType(-3.0, 1.0), seed_gradient(Float64,1,2), seed_gradient(Float64,1,2), false)

   yref_ns_pow1 = MC{2}(16.0, 19.0, EAGO.IntervalType(9.0, 49.0), @SVector[8.0, 0.0], @SVector[10.0, 0.0], false)
   yref_ns_pow2 = MC{2}(0.0625, 0.08843537414965986, EAGO.IntervalType(0.0204081, 0.111112), @SVector[-0.03125, 0.0], @SVector[-0.0226757, 0.0], false)
   yref_ns_pow3 = MC{2}(0.04938271604938271, 0.08246527777777776, EAGO.IntervalType(0.015625, 0.111112), @SVector[0.0219479, 0.0], @SVector[0.0190972, 0.0], false)
   yref_ns_pow4 = MC{2}(20.25, 25.5, EAGO.IntervalType(9.0, 64.0), @SVector[-9.0, 0.0], @SVector[-11.0, 0.0], false)
   yref_ns_pow5 = MC{2}(-172.5, -91.125, EAGO.IntervalType(-512.0, -27.0), @SVector[97.0, 0.0], @SVector[60.75, 0.0], false)
   yref_ns_pow6 = MC{2}(410.0625, 1285.5, EAGO.IntervalType(81.0, 4096.0), @SVector[-364.5, 0.0], @SVector[-803.0, 0.0], false)
   yref_ns_pow7 = MC{2}(91.125, 172.5, EAGO.IntervalType(27.0, 512.0), @SVector[60.75, 0.0], @SVector[97.0, 0.0], false)
   yref_ns_pow8 = MC{2}(410.0625, 1285.5, EAGO.IntervalType(81.0, 4096.0), @SVector[364.5, 0.0], @SVector[803.0, 0.0], false)
   yref_ns_pow9 = MC{2}(410.0625, 2818.5, EAGO.IntervalType(0, 4096.0), @SVector[-364.5, 0.0], @SVector[-365.0, 0.0], false)
   EAGO.set_mc_differentiability!(0); @test check_vs_refv(^, x1, 2, yref_ns_pow1, mctol)
   EAGO.set_mc_differentiability!(0); @test check_vs_refv(^, x1, -2, yref_ns_pow2, mctol)
   EAGO.set_mc_differentiability!(0); @test check_vs_refv(^, x2, -2, yref_ns_pow3, mctol)
   EAGO.set_mc_differentiability!(0); @test check_vs_refv(^, x2, 2, yref_ns_pow4, mctol)
   EAGO.set_mc_differentiability!(0); @test check_vs_refv(^, x2, 3, yref_ns_pow5, mctol)
   EAGO.set_mc_differentiability!(0); @test check_vs_refv(^, x2, 4, yref_ns_pow6, mctol)
   EAGO.set_mc_differentiability!(0); @test check_vs_refv(^, x3, 3, yref_ns_pow7, mctol)
   EAGO.set_mc_differentiability!(0); @test check_vs_refv(^, x3, 4, yref_ns_pow8, mctol)
   EAGO.set_mc_differentiability!(0); @test check_vs_refv(^, x4, 4, yref_ns_pow9, mctol)
   EAGO.set_mc_differentiability!(0); @test x2^(1) == x2

   yref_d1_pow1 = MC{2}(0.25, 0.2857142857142857, EAGO.IntervalType(0.142857, 0.333334), @SVector[-0.0625, 0.0], @SVector[-0.047619, 0.0], false)
   yref_d1_pow2 = MC{2}(16.0, 19.0, EAGO.IntervalType(9.0, 49.0), @SVector[8.0, 0.0], @SVector[10.0, 0.0], false)
   yref_d1_pow3 = MC{2}(16.0, 19.0, EAGO.IntervalType(9.0, 49.0), @SVector[-8.0, 0.0], @SVector[-10.0, 0.0], false)
   yref_d1_pow4 = MC{2}(2.66666666666666, 7.0, EAGO.IntervalType(0.0, 9.0), @SVector[-4.0, 0.0], @SVector[-2.0, 0.0], false)
   yref_d1_pow5 = MC{2}(64.0, 106.0, EAGO.IntervalType(27.0, 343.0), @SVector[48.0, 0.0], @SVector[79.0, 0.0], false)
   yref_d1_pow6 = MC{2}(-106.0, -64.0, EAGO.IntervalType(-343.0, -27.0), @SVector[79.0, 0.0], @SVector[48.0, 0.0], false)
   yref_d1_pow7 = MC{2}(-20.25, -7.750, EAGO.IntervalType(-27.0, 1.0), @SVector[6.75, 0.0], @SVector[12.25, 0.0], false)
   yref_d1_pow8 = MC{2}(0.015625, 0.02850664075153871, EAGO.IntervalType(0.00291545, 0.0370371), @SVector[-0.0117188, 0.0], @SVector[-0.0085304, 0.0], false)
   yref_d1_pow9 = MC{2}(-0.02850664075153871, -0.015625, EAGO.IntervalType(-0.0370371, -0.00291545), @SVector[-0.0085304, 0.0], @SVector[-0.0117188, 0.0], false)
   yref_d1_pow10 = MC{2}(0.00390625, 0.009363382541225106, EAGO.IntervalType(0.000416493, 0.0123457), @SVector[-0.00390625, 0.0], @SVector[-0.0029823, 0.0], false)
   yref_d1_pow11 = MC{2}(0.00390625, 0.009363382541225106, EAGO.IntervalType(0.000416493, 0.0123457), @SVector[0.00390625, 0.0], @SVector[0.0029823, 0.0], false)
   yref_d1_pow12 = MC{2}(256.0, 661.0, EAGO.IntervalType(81.0, 2401.0), @SVector[256.0, 0.0], @SVector[580.0, 0.0], false)
   yref_d1_pow13 = MC{2}(256.0, 661.0, EAGO.IntervalType(81.0, 2401.0), @SVector[-256.0, 0.0], @SVector[-580.0, 0.0], false)
   yref_d1_pow14 = MC{2}(16.0, 61.0,  EAGO.IntervalType(0.0, 81.0), @SVector[-32.0, 0.0], @SVector[-20.0, 0.0], false)
   EAGO.set_mc_differentiability!(1); @test check_vs_refv(^, x1, -1, yref_d1_pow1, mctol)
   EAGO.set_mc_differentiability!(1); @test check_vs_refv(^, x1, 2, yref_d1_pow2, mctol)
   EAGO.set_mc_differentiability!(1); @test check_vs_refv(^, x5, 2, yref_d1_pow3, mctol)
   EAGO.set_mc_differentiability!(1); @test check_vs_refv(^, x6, 2, yref_d1_pow4, mctol)
   EAGO.set_mc_differentiability!(1); @test check_vs_refv(^, x1, 3, yref_d1_pow5, mctol)
   EAGO.set_mc_differentiability!(1); @test check_vs_refv(^, x5, 3, yref_d1_pow6, mctol)
   EAGO.set_mc_differentiability!(1); @test check_vs_refv(^, x6, 3, yref_d1_pow7, mctol)
   EAGO.set_mc_differentiability!(1); @test check_vs_refv(^, x1, -3, yref_d1_pow8, mctol)
   EAGO.set_mc_differentiability!(1); @test check_vs_refv(^, x5, -3, yref_d1_pow9, mctol)
   EAGO.set_mc_differentiability!(1); @test check_vs_refv(^, x1, -4, yref_d1_pow10, mctol)
   EAGO.set_mc_differentiability!(1); @test check_vs_refv(^, x5, -4, yref_d1_pow11, mctol)
   EAGO.set_mc_differentiability!(1); @test check_vs_refv(^, x1, 4, yref_d1_pow12, mctol)
   EAGO.set_mc_differentiability!(1); @test check_vs_refv(^, x5, 4, yref_d1_pow13, mctol)
   EAGO.set_mc_differentiability!(1); @test check_vs_refv(^, x6, 4, yref_d1_pow14, mctol)
end

@testset "Multiplication Operator" begin

    ##### Case 1 #####
    x1 = MC{2}(0.0, 0.0, EAGO.IntervalType(-2.0,1.0), seed_gradient(Float64,1,2), seed_gradient(Float64,1,2), false)
    y1 = MC{2}(1.0, 1.0, EAGO.IntervalType(-1.0,2.0), seed_gradient(Float64,2,2), seed_gradient(Float64,2,2), false)
    yref1 = MC{2}(-1.0, 2.0, EAGO.IntervalType(-4.0,2.0), @SVector[2.0, 1.0], @SVector[2.0, -2.0], false)
    EAGO.set_mc_differentiability!(0); @test check_vs_ref2(*, x1, y1, yref1, mctol)

    ##### Case 2 #####
    x2 = MC{2}(3.0, 3.0, EAGO.IntervalType(1.0,5.0), seed_gradient(Float64,1,2), seed_gradient(Float64,1,2), false)
    y2 = MC{2}(1.0, 1.0, EAGO.IntervalType(-1.0,2.0), seed_gradient(Float64,2,2), seed_gradient(Float64,2,2), false)
    yref2 = MC{2}(1.0, 5.0, EAGO.IntervalType(-5.0,10.0), @SVector[2.0, 5.0], @SVector[2.0, 1.0], false)
    EAGO.set_mc_differentiability!(0); @test check_vs_ref2(*, x2, y2, yref2, mctol)

    ##### Case 3 #####
    x3 = MC{2}(-4.0, -4.0, EAGO.IntervalType(-6.0,-2.0), seed_gradient(Float64,1,2), seed_gradient(Float64,1,2), false)
    y3 = MC{2}(2.0, 2.0, EAGO.IntervalType(1.0,3.0), seed_gradient(Float64,2,2), seed_gradient(Float64,2,2), false)
    yref3 = MC{2}(-10.0, -6.0, EAGO.IntervalType(-18.0,-2.0), @SVector[3.0, -2.0], @SVector[1.0, -2.0], false)
    EAGO.set_mc_differentiability!(0); @test check_vs_ref2(*, x3, y3, yref3, mctol)

    ##### Case 4 #####
    x4 = MC{2}(-4.0, -4.0, EAGO.IntervalType(-6.0,-2.0), seed_gradient(Float64,1,2), seed_gradient(Float64,1,2), false)
    y4 = MC{2}(-5.0, -5.0, EAGO.IntervalType(-7.0,-3.0), seed_gradient(Float64,2,2), seed_gradient(Float64,2,2), false)
    yref4 = MC{2}(16.0, 24.0, EAGO.IntervalType(6.0,42.0), @SVector[-3.0, -2.0], @SVector[-7.0, -2.0], false)
    EAGO.set_mc_differentiability!(0); @test check_vs_ref2(*, x4, y4, yref4, mctol)

    ##### Case 5 #####
    x5 = MC{2}(-4.0, -4.0, EAGO.IntervalType(-6.0,-2.0), seed_gradient(Float64,1,2), seed_gradient(Float64,1,2), false)
    y5 = MC{2}(-5.0, -5.0, EAGO.IntervalType(-7.0,4.0), seed_gradient(Float64,2,2), seed_gradient(Float64,2,2), false)
    yref5 = MC{2}(16.0, 24.0, EAGO.IntervalType(-24.0,42.0), @SVector[-7.0, -6.0], @SVector[-7.0, -2.0], false)
    EAGO.set_mc_differentiability!(0); @test check_vs_ref2(*, x5, y5, yref5, mctol)

    ##### Case 6 #####
    x6 = MC{2}(-2.0, -2.0, EAGO.IntervalType(-3.0,4.0), seed_gradient(Float64,1,2), seed_gradient(Float64,1,2), false)
    y6 = MC{2}(3.0, 3.0, EAGO.IntervalType(1.0,4.0), seed_gradient(Float64,2,2), seed_gradient(Float64,2,2), false)
    yref6 = MC{2}(-8.0, -5.0, EAGO.IntervalType(-12.0, 16.0), @SVector[1.0, -3.0], @SVector[4.0, -3.0], false)
    EAGO.set_mc_differentiability!(0); @test check_vs_ref2(*, x6, y6, yref6, mctol)

    ##### Case 7 #####
    x7 = MC{2}(-2.0, -2.0, EAGO.IntervalType(-3.0,4.0), seed_gradient(Float64,1,2), seed_gradient(Float64,1,2), false)
    y7 = MC{2}(-4.0, -4.0, EAGO.IntervalType(-5.0,-3.0), seed_gradient(Float64,2,2), seed_gradient(Float64,2,2), false)
    yref7 = MC{2}(7.0, 9.0, EAGO.IntervalType(-20.0, 15.0), @SVector[-5.0, -3.0], @SVector[-3.0, -3.0], false)
    EAGO.set_mc_differentiability!(0); @test check_vs_ref2(*, x7, y7, yref7, mctol)

    ##### Testing for Smooth Standard Mult #####
    EAGO.set_mc_differentiability!(1)

    seed1 = seed_gradient(Float64,1,2)
    seed2 = seed_gradient(Float64,2,2)
    x1 = MC{2}(0.0,0.0,EAGO.IntervalType(-200.0,200.0),seed1,seed1,false)
    y1 = MC{2}(200.0,200.0,EAGO.IntervalType(0.0,400.0),seed2,seed2,false)
    z1 = x1*y1
    @test isapprox(z1.cc,40000,atol=1E-4)
    @test isapprox(z1.cv,-40000,atol=1E-4)

    x2 = MC{2}(170.0,170.0,EAGO.IntervalType(100.0,240.0),seed1,seed1,false)
    y2 = MC{2}(250.0,250.0,EAGO.IntervalType(100.0,400.0),seed2,seed2,false)
    z2 = x2*y2
    @test isapprox(z2.cc,53000,atol=1E-4)
    @test isapprox(z2.cv,32000,atol=1E-4)

    x3 = MC{2}(-200.0,-200.0,EAGO.IntervalType(-300.0,-100.0),seed1,seed1,false)
    y3 = MC{2}(-300.0,-300.0,EAGO.IntervalType(-400.0,-200.0),seed2,seed2,false)
    z3 = x3*y3
    @test isapprox(z3.cc,70000,atol=1E-4)
    @test isapprox(z3.cv,50000,atol=1E-4)

    # CHECK ME AGAIN???? -47187.5 new, -47460.9375 old
    x4 = MC{2}(150.0,150.0,EAGO.IntervalType(100.0,200.0),seed1,seed1,false)
    y4 = MC{2}(-250.0,-250.0,EAGO.IntervalType(-500.0,-100.0),seed2,seed2,false)
    z4 = x4*y4
    @test isapprox(z4.cc,-30000,atol=1E-3)
    @test isapprox(z4.cv,-47187.5,atol=1E-3)

    x5 = MC{2}(-150.0,-150.0,EAGO.IntervalType(-200.0,-100.0),seed1,seed1,false)
    y5 = MC{2}(300.0,300.0,EAGO.IntervalType(200.0,400.0),seed2,seed2,false)
    z5 = x5*y5
    @test isapprox(z5.cv,-50000,atol=1E-4)
    @test isapprox(z5.cc,-40000,atol=1E-4)
end

@testset "Division" begin
    EAGO.set_mc_differentiability!(0)
    xIBox = SVector{2,EAGO.IntervalType}([EAGO.IntervalType(-3.0,4.0);EAGO.IntervalType(-5.0,-3.0)])
    a = seed_gradient(Float64,1,2)
    b = seed_gradient(Float64,2,2)
    X = MC{2}(-2.0,-2.0,xIBox[1],a,a,false)
    Y = MC{2}(-4.0,-4.0,xIBox[2],b,b,false)
    out = X/Y
    @test isapprox(out.cc,0.6,atol=1E-6)
    @test isapprox(out.cv,0.41666666,atol=1E-6)
    @test isapprox(out.cc_grad[1],-0.2,atol=1E-6)
    @test isapprox(out.cc_grad[2],0.2,atol=1E-6)
    @test isapprox(out.cv_grad[1],-0.333333,atol=1E-6)
    @test isapprox(out.cv_grad[2],0.1875,atol=1E-6)
    @test isapprox(out.Intv.lo,-1.33333333,atol=1E-6)
    @test isapprox(out.Intv.hi,1.0,atol=1E-6)
end
