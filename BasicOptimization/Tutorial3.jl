#=
Tutorial3.jl
Jacob Child
August 6, 2022
Trying to use the Evolutionary.jl package

=#
using Evolutionary
include("TutorialFunctionFile.jl")
#f(x) = (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2 # Rosenbrock
A = 10.0
n = 2.0
#f(x) = A*n .+ [(x[1])^2 - A*cos(2*pi*x[1])] + [(x[2])^2 - A*cos(2*pi*x[2])] #Rastrigin function
#f(x)=(x[1]+2x[2]-7)^2+(2x[1]+x[2]-5)^2 # Booth
#x0 = [-5.12, 5.12] # for Rosenbrock function initial guess
#x0 = [1., 1.]
#lower = [-5.12, -5.12] # for Rastrigin function lower bound
#upper = [5.12, 5.12] # for Rastrigin function upper bound
lx = [-1.0, 0.0]
ux = [10.0, 10.0]
lc = [-5.0] # lower constraint bound
uc = [Inf]
x0 = [5.0, 2.0]
cnst = WorstFitnessConstraints(lx, ux, lc, uc, ConstraintFunction)
#ga = GA(populationSize=100,selection=uniformranking(3), mutation=gaussian(),crossover=uniformbin())
#opts = Evolutionary.Options(iterations = 1000, abstol = .00000000001, reltol = .00000000001)
result = Evolutionary.optimize(ObjectiveFunction, cnst, GA(populationSize = 1000, mutationRate = .1))
#results = Evolutionary.optimize(ObjectiveFunction, BoxConstraints(lx, ux), GA(populationSize = 1000, mutationRate = .1))