module SNOPT_TestSuite

using Compat
using Compat.Test
using IntervalArithmetic
using EAGO

# PLACE HOLDER FOR SNOPT TESTS

# Tests are set to be formally ignored in coverage since SNOPT isn't open source
# Try running the SNOPT examples to test that a build configured for SNOPT is
# working or uncommenting the below section for local testing

@testset "SNOPT Functionality Testing" begin
y = [1.0,2.0]
YI = [Interval{Float64}(0.0,3.0),Interval{Float64}(0.0,3.0)]
YIMc = [MCInterval{Float64}(0.0,3.0),MCInterval{Float64}(0.0,3.0)]
opt = param = pmid = []
@test_broken snopt_callback_LBD_Imp(y,YI,opt,param,pmid)
@test_broken snopt_callback_LBD_Imp(y,YIMc,opt,param,pmid)
@test_broken SNOPT_LBD_Imp(YI,Int64(2),Int64(3),opt,Float64(1.1))
@test_broken SNOPT_LBD_Imp(YIMc,Int64(2),Int64(3),opt,Float64(1.1))

@test_broken snopt_callback_LBD(y,YI,opt)
@test_broken snopt_callback_LBD(y,YIMc,opt)
@test_broken SNOPT_LBD(YI,Int64(2),Int64(3),opt,Float64(1.1))
@test_broken SNOPT_LBD(YIMc,Int64(2),Int64(3),opt,Float64(1.1))
end

end
