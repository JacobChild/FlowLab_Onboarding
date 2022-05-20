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
using XFoil, Printf, Match

userinput1 = 1; #1 is the default "Match case"
println("Input the 4 NACA numbers ie 2412 \nNACA ")
userinput1 = readline()
userinput1 = parse(Int64, userinput1)

#Beginning of Match "Switch" Statement
@match userinput1 begin
    1 =>
    println("This is the default spot!")
end
