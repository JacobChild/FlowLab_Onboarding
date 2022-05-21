#=
ExtenderNSmoother.jl
Jacob Child
May 20th, 2022

High Level Project Overview-
Files in Project: ExtenderNSmoother.jl, PolarPlotter.jl, AirFoilCoordinatesFunction.jl 
PolarPlotter will be created and made to run as a stand alone file ("default" mode discussed below).
It will then be changed to call AirFoilCoordinatesFunction if need be, and AirFoilCoordinatesFunction
will be turned into a function rather than a stand alone file. ExtenderNSmoother will then be made to
run as it calls PolarPlotter, now turned into a function, for data.

Pseudocode-
Turn Polar Plotter into a function, then using the Dierckx package and the 1-D Spline function 
smooth the already produced data from PolarPlotter. With the smoothed data use the CCBlade Viterna
function to extrapolate the airfoil data from -pi to pi. Check which order of operations works better,
smoothing and then extrapolating or extrapolating then smoothing.

Updated Pseudocode
The Extending function goes from -pi to pi and thus must also take in radians, so convert alpha to rads 

=#

#Packages needed for ALL files in Project to run
using Xfoil, CCBlade, Printf, DataFrames, Plots, FLOWMath, Dierckx
include("PolarPlotter.jl")

#Import Needed Data
alpha, Cl, Cd, Cdp, Cm = PolarPlotter()
alpharads = alpha*pi/180    #alpha needs to be in rads for the viterna extension function

#To see the difference 1/4
ClPlot1 = plot(alpharads, Cl)


#Smoothing 
Spline1 = Spline1D(alpharads, Cl; w=ones(length(alpharads)), k=3, bc="nearest", s=0.017)
SmoothCl = evaluate(Spline1, alpharads)

#To see the difference 2/4
SmoothClPlot = plot(alpharads, SmoothCl)

#Extrapolating/Extending
cr75 = 1    #chord/Rtip at 75%Rtip? Rtip = tip radius
newalpharads, ExtendedCl, ExtendedCd = viterna(vec(alpharads), SmoothCl, Cd, cr75)

#To see the difference 3/4
ExtendedClPlot = plot(newalpharads, ExtendedCl)

#Smoothing Extended Data
Spline2 = Spline1D(newalpharads, ExtendedCl; w=ones(length(newalpharads)), k=3, bc="nearest", s=0.15)
SmoothClExtended = evaluate(Spline2, newalpharads)

#To see the difference 4/4
#ExtendedSmoothClPlot = plot(newalpharads, SmoothClExtended)
ExtendedSmoothClPlot = plot!(newalpharads, SmoothClExtended)
display(ExtendedClPlot)

#Combined Plotter
#plot(ClPlot1, SmoothClPlot, ExtendedClPlot, ExtendedSmoothClPlot, layout = (1,4), legend = false)
