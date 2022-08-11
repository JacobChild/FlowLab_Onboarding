using Evolutionary
# x, y ≥ 0
lx = [0.0, 0.0] # lower bound for values
ux = [Inf, Inf] # upper bound for values
# √xy ≥ 100
c(x) = [ prod(map(e->sqrt(e>0 ? e : 0.0), x)) ] # negative values are zeroed
lc   = [100.0] # lower bound for constraint function
uc   = [ Inf ]   # upper bound for constraint function
con = WorstFitnessConstraints(lx, ux, lc, uc, c)
f(x) = 3x[1]+9x[2] # fitness function
x0 = [1., 1.] # individual
ga = GA(populationSize=100,selection=tournament(7),
        mutation=gaussian(),crossover=intermediate(2))
results = Evolutionary.optimize(f, con, x0, ga)