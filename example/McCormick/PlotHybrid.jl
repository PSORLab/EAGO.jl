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

xl = 1.3
xu = 1.7
intv = 0.01
nx = 41
count2 = 1

tempbox = SVector{1,Interval{Float64}}(1.0*[Interval(xl,xu)])
tempmid =  mid.(tempbox)

Vals = zeros(Float64,nx,1)
x_grid = zeros(Float64,nx)
y_grid = zeros(Float64,nx)
values = zeros(Float64,nx)
ccval = zeros(Float64,nx)
cvval = zeros(Float64,nx)
envccval = zeros(Float64,nx)
envcvval = zeros(Float64,nx)

for i=xl:intv:xu
    temp1 = SMCg{1,Interval{Float64},Float64}(Float64(i),Float64(i),seed_g(Float64,1,1),
                                              seed_g(Float64,1,1),Interval(xl,xu),
                                              false)

    x_grid[count2] = i # y
    values[count2] = (i^2-i)*log(i) # x*(x-1)*log*=(x)
    temp_MC = temp1*(temp1-1)*log(temp1)# relax x*y
    ccval[count2] = temp_MC.cc
    cvval[count2] = temp_MC.cv
    envccval[count2] = (i^2-i)*log(i)
    envcvval[count2] = EAGO.line_seg(i,xl,(xl^2-xl)*log(xl),xu,(xu^2-xu)*log(xu))
    count2 += 1
end

plot(x_grid,[cvval,ccval,values],title = "Standard McCormick Relaxation", label=["Convex McCormick Relaxation" "Concave McCormick Relaxation" "f(x) = (x^2-x)*log(x)"],
xlim=[1.25, 1.75],ylim=[-0.1,0.8],ylabel="f(x)",xlabel="x",lw=3)
gui()

#=
plot(x_grid,[envcvval,envccval,values],title = "Envelope Relaxation", label=["Concave Envelope Relaxation" "Convex Envelope Relaxation" "f(x) = (x^2-x)*log(x)"],
xlim=[1.25, 1.75],ylim=[-0.1,0.8],ylabel="f(x)",xlabel="x",lw=3)
gui()

png("EnvelopePlot")
=#
