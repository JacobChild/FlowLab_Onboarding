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

Sources:
Xfoil.jl Documentation - https://flow.byu.edu/Xfoil.jl/dev/examples/
Flow Math Documentation (ended up not using their smooth function) - http://flow.byu.edu/FLOWMath.jl/stable/#Smoothing-1
Dierckx Documentation (for smoothing) - https://www.juliapackages.com/p/dierckx
CCBlade Documentation (Viterna function used for extending) - https://flow.byu.edu/CCBlade.jl/stable/reference/

    cd("C:/Users/child/Documents/Flow_Lab/Onboarding/XfoilProject")
=#

#Packages needed for ALL files in Project to run
using Xfoil, CCBlade, Printf, DataFrames, Plots, FLOWMath, Dierckx, DelimitedFiles
include("PolarPlotter.jl")

#functions
function smoother(Coef,alpharads, amount)
    #From the Dierckx package, creates a 1-D spline to the kth order with an "s" smoothing amount.
    #bc = boundary conditions and specifies what to do outside of the min/max range, it can be *nearest*, zero, extrapolate, or error
    Spline1 = Spline1D(alpharads, Coef; w=ones(length(alpharads)), k=3, bc="nearest", s=amount) 
    SmoothCoef = evaluate(Spline1, alpharads)
    return SmoothCoef
end

#Import Needed Data
alpha, Cl, Cd, Cdp, Cm = PolarPlotter() #Calls my PolarPlotter function/file which uses Xfoil to get the needed data.
alpharads = alpha.*(pi/180)   #alpha needs to be in rads for the viterna extension function

#To see the difference 1/4
#Parts 1-4 were used during original coding to compare differences in plots
#not needed in final version, but kept for future troubleshooting 
ClPlot1 = plot(alpharads, Cl)
CdPlot1 = plot(alpharads, Cd) 
CdpPlot1= plot(alpharads, Cdp)
CmPlot1 = plot(alpharads, Cm)

 
#Smoothing 
#This calls the smoother function to smooth out the current data, the 3rd input is the smoothing amount
SmoothCl = smoother(Cl,alpharads,.0075)
SmoothCd = smoother(Cd,alpharads,.0065)
SmoothCdp = smoother(Cdp,alpharads,.005)
SmoothCm = smoother(Cm,alpharads,.005)


#To see the difference 2/4
SmoothClPlot = plot(alpharads, SmoothCl)
SmoothCdPlot = plot(alpharads, SmoothCd)
SmoothCdpPlot = plot(alpharads, SmoothCdp)
SmoothCmPlot = plot(alpharads, SmoothCm)

#Extrapolating/Extending
cr75 = 1.0    #chord/Rtip at 75%Rtip? Rtip = tip radius
#Uses the viterna function from the CCBlade package to extend the values we have from -pi to pi
newalpharads, ExtendedCl, ExtendedCd = viterna(vec(alpharads), SmoothCl, SmoothCd, cr75)
#The viterna method was not made for Cdp and Cm, thus those are not Extended

#To see the difference 3/4
ExtendedClPlot = plot(newalpharads, ExtendedCl)

#Smoothing Extended Data
ExtendedSmoothCl = smoother(ExtendedCl,newalpharads,.001)
ExtendedSmoothCd = smoother(ExtendedCd,newalpharads,.005)

#To see the difference 4/4
ExtendedSmoothClPlot = plot!(newalpharads, ExtendedSmoothCl)

#Unit Conversion from radians to degrees, this is optional and could be dangerous as I am not changing the variable name
alpharads = alpharads*180/pi
#newalpharads = newalpharads*180/pi

"""
Matrix Meddling
This is to combine and arrange the data to prepare for plotting
"""
#Combines the alpharads (original and extended) arrays into 2 column vectors
Combalpharads = cat(alpharads,alpharads, dims = (2,2))  #I don't fully understand what the dims does, but it works that way
Combnewalpharads = cat(newalpharads, newalpharads, dims = (2,2))

#combines the smoothed and non smoothed Cl values into 2 Col vectors 
CombCl = cat(Cl, SmoothCl, dims = (2,2))
CombExtendedCl = cat(ExtendedCl, ExtendedSmoothCl, dims = (2,2)) 

# as above, but Cd
CombCd = cat(Cd, SmoothCd, dims = (2,2))
CombExtendedCd = cat(ExtendedCd, ExtendedSmoothCd, dims = (2,2)) 

#As above, but Cdp
CombCdp = cat(Cdp, SmoothCdp, dims = (2,2))

#As above, but Cm
CombCm = cat(Cm, SmoothCm, dims = (2,2))


"""
Plot Palooza
Saves the various plots as plot objects and combines all plot objects of a coefficient into one plot object
"""
CombClPlot = plot(Combalpharads,CombCl, title = "Cl Plots", label = ["Org" "Smooth"])
CombExtendedClPlot = plot(Combnewalpharads, CombExtendedCl, label = ["Org" "Smooth"])
FullClPlot = plot(CombClPlot,CombExtendedClPlot, layout = (2,1), legend = false)

CombCdPlot = plot(Combalpharads,CombCd, title = "Cd Plots", label = ["Org" "Smooth"])
CombExtendedCdPlot = plot(Combnewalpharads, CombExtendedCd, label = ["Org" "Smooth"])
FullCdPlot = plot(CombCdPlot,CombExtendedCdPlot, layout = (2,1), legend = false)

CombCdpPlot = plot(Combalpharads,CombCdp, title = "Cdp Plots", label = ["Org" "Smooth"])

CombCmPlot = plot(Combalpharads,CombCm, title = "Cm Plots", label = ["Org" "Smooth"])

#Creates an array of Plot objects (from the above Plot Palooza) and plots them
PlotArray = []
push!(PlotArray,FullClPlot)
push!(PlotArray,FullCdPlot)
push!(PlotArray,CombCdpPlot)
push!(PlotArray,CombCmPlot)
plot(PlotArray...)
plot!(size = (1000,900),legend = true)    #changes the plot size 

#Output to a file
#touch("EpplerE63Data.txt")  #creates the file
open("EpplerE63Data.txt", "w") do FileID
lines = ("EppelerE63 Data\n")
write(FileID, lines)
lines = ("1000000\n")
write(FileID, lines)
lines = ("0\n")
write(FileID, lines)

#Write table/data to file 
DataTable = cat(newalpharads, ExtendedSmoothCl, ExtendedSmoothCd, dims = (2,2))
writedlm(FileID, DataTable)     #uses DelimitedFiles package, I believe default dlm is \t

end