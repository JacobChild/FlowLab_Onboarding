#=
Tutorial2.jl
Jacob Child
August 6, 2022
Trying to use the Optim.jl package
=#
using Optim, OptimTestProblems

#Setup
#f(x) = (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2 #Rosenbrock function
#A = 10
#n = 2
#f(x) = A*n +[(x[1])^2 - A*cos(2*pi*x[1])] + [(x[2])^2 - A*cos(2*pi*x[2])] #Rastrigin function
#x0 = [0.0, 0.0] #for Rosenbrock function initial guess

prob = UnconstrainedProblems.examples["Rosenbrock"];
#Solve
res = optimize(prob.f, fill(-100.0, 2), fill(100.0, 2), prob.initial_x, SAMIN(), Optim.Options(iterations=10^6))