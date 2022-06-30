#=
TSExtenderNSmootherFunction.jl
Jacob Child
June 29th, 2022

High Level Project Overview-
TradeStudyTopFile.jl will extract needed airfoil data from the given csv files in the data folder.
That data will then be fed into this file (TSExtenderNSmootherFunction) which will call the 
TSPolarPlotterFunction and TSAirFoilCoordinatesFunction to generate needed airfoil data from xfoil.

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

=#
function TSExtenderNSmootherFunction(alpharads, Cl, Cd, Cdp, Cm, cr75)
#functions
function smoother(Coef,alpharads, amount)
    #From the Dierckx package, creates a 1-D spline to the kth order with an "s" smoothing amount.
    #bc = boundary conditions and specifies what to do outside of the min/max range, it can be *nearest*, zero, extrapolate, or error
    Spline1 = Spline1D(alpharads, Coef; w=ones(length(alpharads)), k=3, bc="nearest", s=amount) 
    SmoothCoef = evaluate(Spline1, alpharads)
    return SmoothCoef
end
 
#Smoothing 
#This calls the smoother function to smooth out the current data, the 3rd input is the smoothing amount
SmoothCl = smoother(Cl,alpharads,.0075)
SmoothCd = smoother(Cd,alpharads,.0065)
SmoothCdp = smoother(Cdp,alpharads,.005)
SmoothCm = smoother(Cm,alpharads,.005)

#Extrapolating/Extending
#Uses the viterna function from the CCBlade package to extend the values we have from -pi to pi
newalpharads, ExtendedCl, ExtendedCd = viterna(vec(alpharads), SmoothCl, SmoothCd, cr75)
#The viterna method was not made for Cdp and Cm, thus those are not Extended

#Smoothing Extended Data
ExtendedSmoothCl = smoother(ExtendedCl,newalpharads,.001)
ExtendedSmoothCd = smoother(ExtendedCd,newalpharads,.005)

#Unit Conversion from radians to degrees
alphadegs = alpharads*180/pi
newalphadegs = newalpharads*180/pi

"""
Matrix Meddling
This is to combine and arrange the data to prepare for plotting
"""
#Combines the alpharads (original and extended) arrays into 2 column vectors
Combalphadegs = cat(alphadegs,alphadegs, dims = (2,2))  #I don't fully understand what the dims does, but it works that way
Combnewalphadegs = cat(newalphadegs, newalphadegs, dims = (2,2))

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
CombClPlot = plot(Combalphadegs,CombCl, title = "Cl Plots", label = ["Org" "Smooth"])
CombExtendedClPlot = plot(Combnewalphadegs, CombExtendedCl, label = ["Org" "Smooth"])
FullClPlot = plot(CombClPlot,CombExtendedClPlot, layout = (2,1), legend = false)

CombCdPlot = plot(Combalphadegs,CombCd, title = "Cd Plots", label = ["Org" "Smooth"])
CombExtendedCdPlot = plot(Combnewalphadegs, CombExtendedCd, label = ["Org" "Smooth"])
FullCdPlot = plot(CombCdPlot,CombExtendedCdPlot, layout = (2,1), legend = false)

CombCdpPlot = plot(Combalphadegs,CombCdp, title = "Cdp Plots", label = ["Org" "Smooth"])

CombCmPlot = plot(Combalphadegs,CombCm, title = "Cm Plots", label = ["Org" "Smooth"])

#Creates an array of Plot objects (from the above Plot Palooza) and plots them
PlotArray = []
push!(PlotArray,FullClPlot)
push!(PlotArray,FullCdPlot)
push!(PlotArray,CombCdpPlot)
push!(PlotArray,CombCmPlot)
plot(PlotArray...)
display(FullClPlot)
#plot!(size = (1000,900),legend = true)    #changes the plot size 

#Output to a file
#touch("CurrentAirfoilData.txt")  #creates the file
open("CurrentAirfoilData.txt", "w") do FileID
lines = ("Current Airfoil Data\n")
write(FileID, lines)
lines = ("$Re\n")
write(FileID, lines)
lines = ("$M\n")
write(FileID, lines)

#Write table/data to file 
DataTable = cat(newalpharads, ExtendedSmoothCl, ExtendedSmoothCd, dims = (2,2)) #remember this needs rads
writedlm(FileID, DataTable)     #uses DelimitedFiles package, I believe default dlm is \t

end

return PlotArray, newalphadegs, ExtendedSmoothCl, ExtendedSmoothCd

end