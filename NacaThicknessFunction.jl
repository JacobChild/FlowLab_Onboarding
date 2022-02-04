#=
NacaThicknessFunction.jl
Jacob Child
February 3, 2022
Description:
NACA 4-series Thickness Formula The NACA 4-series formula for airfoil thickness
along the chord is given by the following polynomial function, where m is the value of maximum
thickness (as a percentage of the chord), and x is the x-position along the chord. (This formulation
is for airfoils with sharp trailing edges.)
t = 10m(.2969sqrt(x) - .1620x - .3537x^2 + .2843x^3 - .1015x^4)

Pseudo code- Have a function to ask for and generate x values in an array
Then have a function to calculate various thicknesses
Format as table? and output to screen
***use docstrings***
=#

"""
xinput(x, step)
Takes in the desired length, l, and step size and outputs an x array

"""
function x_input(L, step)
    x = [0:step:L;] #makes array stopping at L with step size step
    return x
end

"""
thickness_calculations(x)
Takes in the x coordinates, calculates thickness, and ouputs the array

"""
function thickness_calculations(x)
    m = 0.10
    t = 10*m.*(.2969*sqrt.(x) - .1260*x - .3537*x.^2 + .2843*x.^3 - .1015*x.^4)
    return t
end

#using Pkg  #Just for the first run through
#Pkg.add("DataFrames")  #just for the first run through

println("Hello! We are going to calculate the thickness along an airfoil \n")
println("\n Enter the Chord length")
L = readline()
L = parse(Int64, L)

println("\n Enter the desired step size to calculate the thickness at")
step = readline()
step = parse(Float64, step)

#println("Type of Inputs", typeof(L), typeof(step))

x = x_input(L, step)
t = thickness_calculations(x)

#println(t) #If Running in VS Code Uncomment this as I couldn't get DataFrame to work

using DataFrames #Only seems to work in Julia Terminal, not VS Code Terminal? ***research

DataFrame()

DataFrame("X" => x, "Thickness" => t)


