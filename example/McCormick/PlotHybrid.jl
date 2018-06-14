#workspace()
using EAGO
using DataFrames
using IntervalArithmetic
using StaticArrays
#using DataFrames
using CSV
using Plots
#using PyPlot
pyplot()

xl = 1.0
xu = 4.0
nx = 41
intv = (xu-xl)/(nx-1)
count2 = 1

set_diff_relax(1)

tempbox = SVector{1,Interval{Float64}}(1.0*[Interval(xl,xu)])
tempmid =  mid.(tempbox)

Vals = zeros(Float64,nx,1)
x_grid = zeros(Float64,nx)
y_grid = zeros(Float64,nx)
values = zeros(Float64,nx)
ccval = zeros(Float64,nx)
cvval = zeros(Float64,nx)
xLval = zeros(Float64,nx)
xUval = zeros(Float64,nx)
envccval = zeros(Float64,nx)
envcvval = zeros(Float64,nx)

temp_MC = SMCg{1,Interval{Float64},Float64}(Float64(2.5),Float64(2.5),seed_g(Float64,1,1),
                                          seed_g(Float64,1,1),Interval(xl,xu),
                                          false)
temp1 = temp_MC
temp_MC = temp1*(temp1-5.0)*sin(temp1)
refccgrad = temp_MC.cc_grad[1]
refcvgrad = temp_MC.cv_grad[1]
refcc = temp_MC.cc
refcv = temp_MC.cv
for i=xl:intv:xu
    temp1 = SMCg{1,Interval{Float64},Float64}(Float64(i),Float64(i),seed_g(Float64,1,1),
                                              seed_g(Float64,1,1),Interval(xl,xu),
                                              false)

    x_grid[count2] = i # y
    values[count2] = i*(i-5.0)*sin(i) # x*(x-1)*log*=(x)
    temp_MC = temp1*(temp1-5.0)*sin(temp1)# relax x*y
    ccval[count2] = temp_MC.cc
    cvval[count2] = temp_MC.cv
    xLval[count2] = temp_MC.Intv.lo
    xUval[count2] = temp_MC.Intv.hi
    envccval[count2] = refcc + refccgrad*(2.5-i)*0.6
    envcvval[count2] = refcv + refcvgrad*(2.5-i)
    count2 += 1
end
#=
plot(x_grid,[values,cvval,ccval,envcvval,envccval,xLval,xUval],title = "Standard McCormick Relaxation",
label=["f(x) = x*(x-5)*sin(x)" "Convex McCormick Relaxation" "Concave McCormick Relaxation" "Concave Affine Relaxation" "Convex Affine Relaxation" "Convex Natural Interval" "Concave Natural Interval"],
xlim=[1.00, 4.00],ylim=[-20,20],ylabel="f(x)",xlabel="x",lw=3)
gui()
=#

plot(x_grid,[values,cvval,ccval,xLval,xUval],title = "Differentiable McCormick Relaxation",
label=["f(x) = x*(x-5)*sin(x)" "Convex McCormick Relaxation" "Concave McCormick Relaxation" "Convex Natural Interval" "Concave Natural Interval"],
xlim=[1.00, 4.00],ylim=[-20,20],ylabel="f(x)",xlabel="x",lw=3)
gui()
#=
plot(x_grid,[envcvval,envccval,values],title = "Envelope Relaxation", label=["Concave Envelope Relaxation" "Convex Envelope Relaxation" "f(x) = (x^2-x)*log(x)"],
xlim=[1.25, 1.75],ylim=[-0.1,0.8],ylabel="f(x)",xlabel="x",lw=3)
gui()

png("EnvelopePlot")
=#
