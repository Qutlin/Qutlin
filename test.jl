#=
test:
- Julia version: 
- Author: fehse
- Date: 2021-03-26
=#
using Plots
x = 1:10; y = rand(10); # These are the plotting data
plt = plot(x, y)
display(plt)
