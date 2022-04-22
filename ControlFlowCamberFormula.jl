#=
ControlFlowCamberFormula.jl
Jacob Child
Flow Lab Onboarding
Feb 3, 2022
Background: A Camber line is "a measure of the curvature of an airfoil"
    It is halfway between the upper and lower curves of the airfoil

    Equation
    if X<=p
        zbar = c*(2*p*x-x^2)/(p^2)
    else if x > p 
        zbar = c*(1-2*p+2*p*x-x^2)/((1-p)^2)
    else 
        println("error")

Pseudocode - ask for the c and p inputs , run an if statement and have A
for loop run through x values, 0:.1:1;
=#

println("Hello! We are going to calculate the the camber of an airfoil \n")
println("\n input c \n")
input1 = readline()
c = parse(Float64, input1)
println("\n input p \n")
input2 = readline()
p = parse(Float64, input2)
#println(p) #to check if it parsed correctly
x = [0:.1:1.1;]    #makes an array stopping at 1 and step size .1
global zbar = x[1:11]
global k = 1
while x[k] <= 1
    if x[k] <= p
        global zbar[k] = c*(2*p*x[k]-x[k]^2)/(p^2)
        global k = k + 1 #iterate

    elseif x[k] > p 
        global zbar[k] = c*(1-2*p+2*p*x[k]-x[k]^2)/((1-p)^2)
        global k = k + 1 #iterate

    else 
        println("error") 
        break
    end

end
println("\n")
println(zbar)

pop!(x) #removes the last value of the x array to allow DataFrames to work
using DataFrames #only works in Julia REPL?
DataFrame()
DataFrame("X" => x, "zbar" => zbar)