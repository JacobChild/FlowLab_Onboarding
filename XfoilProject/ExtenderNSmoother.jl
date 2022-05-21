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

#functions
function smoother(Coef,alpharads, amount)
    Spline1 = Spline1D(alpharads, Coef; w=ones(length(alpharads)), k=3, bc="nearest", s=amount)
    SmoothCoef = evaluate(Spline1, alpharads)
    return SmoothCoef
end

#Import Needed Data
alpha, Cl, Cd, Cdp, Cm = PolarPlotter()
alpharads = alpha*pi/180    #alpha needs to be in rads for the viterna extension function

#To see the difference 1/4
ClPlot1 = plot(alpharads, Cl)
CdPlot1 = plot(alpharads, Cd) 
CdpPlot1= plot(alpharads, Cdp)
CmPlot1 = plot(alpharads, Cm)


#Smoothing 
SmoothCl = smoother(Cl,alpharads,.017)
SmoothCd = smoother(Cd,alpharads,.017)
SmoothCdp = smoother(Cdp,alpharads,.017)
SmoothCm = smoother(Cm,alpharads,.017)


#To see the difference 2/4
SmoothClPlot = plot(alpharads, SmoothCl)
SmoothCdPlot = plot(alpharads, SmoothCd)
SmoothCdpPlot = plot(alpharads, SmoothCdp)
SmoothCmPlot = plot(alpharads, SmoothCm)

#Extrapolating/Extending
cr75 = 1    #chord/Rtip at 75%Rtip? Rtip = tip radius
newalpharads, ExtendedCl, ExtendedCd = viterna(vec(alpharads), SmoothCl, SmoothCd, cr75)
newalpharads, ExtendedCdp, ExtendedCm = viterna(vec(alpharads), SmoothCdp, SmoothCm, cr75)

#To see the difference 3/4
ExtendedClPlot = plot(newalpharads, ExtendedCl)

#Smoothing Extended Data
ExtendedSmoothCl = smoother(ExtendedCl,newalpharads,.025)
ExtendedSmoothCd = smoother(ExtendedCd,newalpharads,.025)
ExtendedSmoothCdp = smoother(ExtendedCdp,newalpharads,.05)
ExtendedSmoothCm = smoother(ExtendedCm,newalpharads,.05)

#To see the difference 4/4
#ExtendedSmoothClPlot = plot(newalpharads, SmoothClExtended)
ExtendedSmoothClPlot = plot!(newalpharads, ExtendedSmoothCl)
display(ExtendedClPlot)

#Matrix Meddling
Combalpharads = cat(alpharads,alpharads, dims = (2,2))  #I don't fully understand what the dims does, but it works that way
Combnewalpharads = cat(newalpharads, newalpharads, dims = (2,2))

CombCl = cat(Cl, SmoothCl, dims = (2,2))
CombExtendedCl = cat(ExtendedCl, ExtendedSmoothCl, dims = (2,2)) 

CombCd = cat(Cl, SmoothCd, dims = (2,2))
CombExtendedCd = cat(ExtendedCd, ExtendedSmoothCd, dims = (2,2)) 

CombCdp = cat(Cl, SmoothCl, dims = (2,2))
CombExtendedCdp = cat(ExtendedCdp, ExtendedSmoothCdp, dims = (2,2))

CombCm = cat(Cl, SmoothCl, dims = (2,2))
CombExtendedCm = cat(ExtendedCm, ExtendedSmoothCm, dims = (2,2)) 

#Plot Palooza
CombClPlot = plot(Combalpharads,CombCl, title = "Cl Plots")
CombExtendedClPlot = plot(Combnewalpharads, CombExtendedCl)
FullClPlot = plot(CombClPlot,CombExtendedClPlot, layout = (2,1), legend = false)

CombCdPlot = plot(Combalpharads,CombCd, title = "Cd Plots")
CombExtendedCdPlot = plot(Combnewalpharads, CombExtendedCd)
FullCdPlot = plot(CombCdPlot,CombExtendedCdPlot, layout = (2,1), legend = false)

CombCdpPlot = plot(Combalpharads,CombCdp, title = "Cdp Plots")
CombExtendedCdpPlot = plot(Combnewalpharads, CombExtendedCdp)
FullCdpPlot = plot(CombCdpPlot,CombExtendedCdpPlot, layout = (2,1), legend = false)

CombCmPlot = plot(Combalpharads,CombCm, title = "Cm Plots")
CombExtendedCmPlot = plot(Combnewalpharads, CombExtendedCm)
FullCmPlot = plot(CombCmPlot,CombExtendedCmPlot, layout = (2,1), legend = false)

PlotArray = []
push!(PlotArray,FullClPlot)
push!(PlotArray,FullCdPlot)
push!(PlotArray,FullCdpPlot)
push!(PlotArray,FullCmPlot)
plot(PlotArray...)
