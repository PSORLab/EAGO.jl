#workspace()

using EAGO
using DataFrames
using IntervalArithmetic
using StaticArrays
#using DataFrames
using CSV
using Plots
pyplot()
#using PyPlot

EAGO.set_diff_relax(1)

xl = -4.0
xu = 5.0

yl = -5.00
yu = 3.0

intv = 0.2
nx = Integer((xu-xl)/intv+1)
ny = Integer((yu-yl)/intv+1)

x_grid = zeros(Float64,nx,ny)
y_grid = zeros(Float64,nx,ny)
ccval = zeros(Float64,nx,ny)
cvval = zeros(Float64,nx,ny)
val = zeros(Float64,nx,ny)

count1 = 1
for j=xl:intv:xu
    temp1 = SMCg{2,Interval{Float64},Float64}(Float64(j),Float64(j),seed_g(Float64,1,2),
                                                  seed_g(Float64,1,2),Interval(xl,xu),
                                                  false)
    count2 = 1
    for i=yl:intv:yu
        temp2 = SMCg{2,Interval{Float64},Float64}(Float64(i),Float64(i),seed_g(Float64,1,2),
                                                      seed_g(Float64,1,2),Interval(yl,yu),
                                                      false)
        x_grid[count1,count2] = j # x
        y_grid[count1,count2] = i # y
        outer = max(temp1,temp2)
        ccval[count1,count2] = outer.cc
        cvval[count1,count2] = outer.cv
        val[count1,count2] = max(j,i)
        count2 +=1
    end
    count1 += 1
end


surface(x_grid,y_grid,cvval,xlab = "x", ylab = "y", zlab = "cv")
surface!(x_grid,y_grid,ccval, xlab = "x", ylab = "y", zlab = "cc")
surface!(x_grid,y_grid,val, xlab = "x", ylab = "y", zlab = "val",
         title = "Differentiable McCormick Relaxation")
gui()
