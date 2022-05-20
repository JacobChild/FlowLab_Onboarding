#=
PolarPlotter.jl
Jacob Child
May 20th, 2022

High Level Project Overview-
Files in Project: ExtenderNSmoother.jl, PolarPlotter.jl, AirFoilCoordinatesFunction.jl 
PolarPlotter will be created and made to run as a stand alone file ("default" mode discussed below).
It will then be changed to call AirFoilCoordinatesFunction if need be, and AirFoilCoordinatesFunction
will be turned into a function rather than a stand alone file. ExtenderNSmoother will then be made to
run as it calls PolarPlotter, now turned into a function, for data.

Pseudocode-
Purpose- This code will generate the polar plots of the given Airfoil. It will ultimatley
probably export those results to the ExtenderNSmoother

Layout-
The user will choose whether to run from a loaded coordinates file, plot a new one, or use the default.
This can be done through state code with 3 states. Once the coordinates have been imported, XFoil will
be called and used  to generate the polar plots from -10deg to 15deg.

=#

#Needed Packages
using Xfoil, Printf, DataFrames, Plots

#FUNCTIONS
function readfile(filename)
    x = Float64[]
    y = Float64[]
    open(filename, "r") do file   #opens the file in read mode "r", "do f" creates a do block so the file is closed at the end 
        for line in Iterators.drop(eachline(file), 1)   #Iterators.drop, drops the first iterator as that is the title line
            CurrentLine = split(chomp(line)) #chomp gets rid of the "newline" (might not be necessary),split splits the line with space as the delimiter
            push!(x, parse(Float64,CurrentLine[1]))
            push!(y, parse(Float64, CurrentLine[2]))
        end
    end
    return x, y
end

function fluidconditions()
    #Uncomment the comment block to ask for user input
    #And delete the below
    Re = 10000
    mach = 0
    ncrit = 7
    alpha = -10:1:15
    #=
    println("Input the Reynolds Number\n")
    Re = readline()
    Re = parse(Float64,Re)
    println("Input the Mach Number\n")
    mach = readline()
    mach = parse(Float64, mach)
    println("Input the Ncrit multiplier")
    ncrit = readline()
    ncrit = parse(Float64, ncrit)
    println("Assuming -10 to 15 deg with a 1 deg step size")
    alpha = -10:1:15    #from -10 to 15 with a step size of 1
    =#
    return Re, mach, alpha, ncrit

end




#MAIN CODE
userinput1 = 1; #1 is the default "Match case"
#Uncomment to run other than default

println("Input one of the following\n1 Default NACA 2412 \n2 Load File \n3 Create Profile ")
userinput1 = readline()
userinput1 = parse(Int64, userinput1)


#Beginning of if statement for "mode" of running
x = 0
y = 0
if userinput1 == 1 
    filename = "NACA2412.dat"
    x, y = readfile(filename)

elseif userinput1 == 2
    println("Enter the Filename \n")
    filename = readline()
    x, y = readfile(filename)

elseif userinput1 == 3
    include("AirFoilCoordinatesFunction.jl")
    x, y = AirFoilCoordinatesFunction()

else
    println("Something went wrong")

end

#Solve for the Polars
Xfoil.set_coordinates(x,y)
Xfoil.pane()    #this "panes" the coordinates, not sure what it does, smooths I think?
Re, mach, alpha, ncrit = fluidconditions()

#OUTPUTS
n = length(alpha)
Cl = zeros(n)
Cd = zeros(n)
Cdp = zeros(n) 
Cm = zeros(n)
Converged = zeros(Bool, n)

#Calculations 
for i = 1:n
    Cl[i], Cd[i], Cdp[i], Cm[i], Converged[i] = Xfoil.solve_alpha(alpha[i], Re; mach, iter = 500, ncrit)
    #The above calculates the coefficients at each angle (i of alpha)
    #for some reason it isn't converging on the last one?
end

#Present Data
DataFrame()

Data = DataFrame("Angle" => alpha, "Cl"=> Cl, "Cd" => Cd, "Cdp" => Cdp, "Cm" => Cm, "Convergence" => Converged)

display(Data)

#Couldn't get the plot to work?
Plot1 = plot(x,y)
plot!(xlims = (-.1, 1.1), title = ("NACA Airfoil"), xlabel = "Position", ylabel = "Position")
plot!(aspectratio = :equal)
println(x)
display(Plot1)