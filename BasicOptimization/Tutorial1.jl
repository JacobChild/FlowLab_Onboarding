#=
Tutorial1.jl
Jacob Child
August 6, 2022

Following the tutorial found https://jump.dev/JuMP.jl/stable/tutorials/getting_started/getting_started_with_JuMP/

=#
using JuMP, HiGHS, Ipopt

#Setup 
#model = Model(HiGHS.Optimizer)
model = Model(Ipopt.Optimizer)
@variable(model, -5.12 <= x <= 5.12)
@variable(model, -5.12 <= y <= 5.12)
A = 10
n = 2
f(x,y) = A*n +[(x)^2 - A*cos(2*pi*x)] + [(y)^2 - A*cos(2*pi*y)] #Rastrigin function
register(model, f, 2, f)
@NLconstraint(model, -5.12 <= x <= 5.12)
@NLconstraint(model, -5.12 <= y <= 5.12)

@NLobjective(model, Min, f(x,y)) #Rastrigin function

#See where we are at
print(model)

#Solve the model
optimize!(model)

#! I am bailing on using JuMP because of the following from their documentation
#=
JuMP does support nonlinear programs with constraints and objectives containing user-defined functions.
However, the functions must be automatically differentiable, or need to provide explicit derivatives.
(See User-defined Functions for more information.)
If your function is a black-box that is non-differentiable (e.g., the output of a simulation written
in C++), JuMP is not the right tool for the job. This also applies if you want to use a derivative
free method.
Even if your problem is differentiable, if it is unconstrained there is limited benefit
(and downsides in the form of more overhead) to using JuMP over tools which are only concerned with
function minimization.

Alternatives to consider are:

Optim.jl
GalacticOptim.jl
NLopt.jl
=#