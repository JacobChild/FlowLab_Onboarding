#TradeStudyTopFile.jl
#Jacob Child
#June 28th, 2022
#= Purpose: To compare the Coefficients of Thrust, Torque, and Power when pitch, chord, 
advance ratio, and radius are independently varied.

Pseudocode: Open the high level csv file and from that get the neccessary propeller data 
and file names for where to extract the addtional data.

 cd("C:/Users/child/Documents/Flow_Lab/Onboarding/TradeStudy")

 Note: I used the "Better Comments" Extension to help find things, the legend is below
 * for highlighted text
! for errors and warnings
? for queries and questions
// for strikethrough
TODO for to-dos
=#

#Needed packages for the whole project, and function file 
using Xfoil, CCBlade, DataFrames, Plots, FLOWMath, Dierckx, DelimitedFiles, Plots.Plots.Measures
#Plotting Settings
dims = (1,2) #for plot layouts after each parameter that was varied, ie "subplot" layouts
pad = 16mm #padding/margin for the plots
#pyplot() #plot package that is being used so the plots come up as seperate figures in the VSCode plot navigator pane 
using PyCall #allows figures you can interact with, ie rotate
#pygui(true) #makes it so I can rotate the plots
using PyPlot
Plots.default(reuse = false) #means that each figure (in pyplot) doesn't overwrite the one before it

#File Parameters
NumofVars = 50 #length of the variable vectors, ie if =20 then 20 pitches will be tested

#Import the high level csv file
CmdFileName = "DJI-II.csv"
cmdfile = readdlm("data\\rotors\\$CmdFileName", ',')
Rtip = cmdfile[2,2]
#Rhub = 0.004938276 #This is the first r value in the DJI-II prop, there is an interference, but that is fixed below
Rhub = cmdfile[3,2]
NumofBlades = cmdfile[4,2]
BladeFileName = cmdfile[5,2]

#Define the rotor object
rotor = Rotor(Rhub, Rtip, NumofBlades) #stores a "Rotor" object

#Import the blade directory file
BladeFile = readdlm("data\\rotors\\$BladeFileName", ',')

ChordFile = BladeFile[2,2]
PitchFile = BladeFile[3,2]
SweepFile = BladeFile[4,2]

HeightFile = BladeFile[5,2]
AirfoilFiles = BladeFile[6,2]
SplineOrder = BladeFile[7,2] #?Is that what was used to smooth the data, or what I should use to smooth the data?
SplineSmoothing = BladeFile[8,2] #? Is that was was used to smooth the data, or what I should use to smooth the data?

#Import the ChordFile 
ChordDist = readdlm("data\\rotors\\$ChordFile", ',')
r = ChordDist[2:end,1] .* Rtip 
chord = ChordDist[2:end,2] .* Rtip

#Check if the Geometry is okay
#! if the first value of r < Rhub it will throw the acos bounds error
if r[1] < Rhub
    Rhub = r[1]
    rotor = Rotor(Rhub, Rtip, NumofBlades) #redefines and stores the "Rotor" object
    println("RHub was redefined to: $Rhub, and the rotor object was redefined")
end

#Import the PitchFile
PitchDist = readdlm("data\\rotors\\$PitchFile", ',')
twist = PitchDist[2:end,2] .* pi/180    #import twist and convert to radians

#Import the AirfoilFiles 
AirfoilFiles = readdlm("data\\rotors\\$AirfoilFiles", ',')
Locations = AirfoilFiles[2:end,1] .* Rtip
NumofSections = length(Locations)
ContourFiles = AirfoilFiles[2:end,2] #These are the files that contain the airfoil coordinates
AeroFiles = AirfoilFiles[2:end,3]  #These files have the xfoil output etc, but aren't extrapolated
#AeroFiles may not exist depending on the csv file and whether or not the data was previously calculated
#If above is the case, check airfoiltools.com, or generate your own data with ExtenderNSmoother.jl 
#     #found in the XfoilProject folder 

#Import the ContourFiles
ContourX = [Vector{Float64}(undef,(0)) for _ in 1:length(ContourFiles)] #Initialize the array of vectors
ContourY = [Vector{Float64}(undef,(0)) for _ in 1:length(ContourFiles)] #Initialize the array of vectors

#For loop that extracts the x and y coordinates from each of the ContourFiles
for i = 1:NumofSections
    GeomFileName = ContourFiles[i]
    GeomFiles = readdlm("data\\airfoils\\$GeomFileName", ',')
    ContourX[i] = GeomFiles[2:end,1]
    ContourY[i] = GeomFiles[2:end,2]
end

#Define the fluid conditions and other values for functions
#cr75
indexval = round(Int64, .75*length(chord))
cr75 = chord[indexval]/Rtip

#Fluid Conditions 
rho = 1.225 #kg/m^3
mu = 1.78E-5 #kg/m/s
Vinf = 0  #inflow velocity, this means that the solve done below is for a singular advance ratio (hover in this case)
RPM = 5000
Omega = RPM*pi/30  #converts rpm to rad/s, derivation: rpm*pi*360deg/(60sec*180deg)-> rpm*pi/30


#Includes needed files for the functions
include("TSPolarPlotterFunction.jl") 
include("TSExtenderNSmootherFunction.jl")
#Prep for Loop
aftypes = Array{AlphaAF}(undef, NumofSections) #Initialize the array of alphaAF objects

#For Loop that extracts and uses the data from the aerofiles to create the airfoil objects
for i in 1:NumofSections
    CurrentData = readdlm("data\\airfoils\\$(AeroFiles[i])", ',')
    CurrentHeaderData = CurrentData[1:9,:]
    CurrentAeroData = Float64.(CurrentData[12:end,:])
    Re = CurrentHeaderData[4,2]
    M = CurrentHeaderData[6,2]
    Ncrit = CurrentHeaderData[5,2]
    alphadegs = CurrentAeroData[:,1]
    alpharads = alphadegs * pi/180
    Cl = CurrentAeroData[:,2]
    Cd = CurrentAeroData[:,3]
    Cdp = CurrentAeroData[:,4]
    Cm = CurrentAeroData[:,5]
    OldPlotArray, newalphadegs, ExtendedSmoothCl, ExtendedSmoothCd = TSExtenderNSmootherFunction(alpharads, Cl, Cd, Cdp, Cm, cr75, Re, M)
    aftypes[i] = AlphaAF("CurrentAirfoilData.txt");   #this gives angle of attack, lift coefficient, drag coefficient
end


#For loop to put which airfoil is used for each section into the af_idx array
af_idx = zero(r) #? Why does zero(r) seem to work in this case, but EffRVar = zero(NumofRadiis) doesn't work
#? Ans: Its because zero is based off of length and makes an array, so one input = one zero, unless that 
#?     input is a vector, then it makes that many zeros, it can also take in array dimensions ie zero([3,2])
#?     zeros is specifically to create vectors of a certain length and type, and cannot do arrays (that I know of)

for k = 1:NumofSections
	for i = 1:length(r) 
		if (k != NumofSections)
			if (r[i] >= Locations[k] && r[i] < Locations[k+1])
				af_idx[i] = k

			elseif (r[i] >= Locations[k] && r[i] > Locations[k+1])
				break
			
			end
		elseif(k == NumofSections && r[i] >= Locations[k])
			af_idx[i] = k
		end
	end
end

af_idx = Int.(af_idx)   #Converts the index values to Ints
airfoils = aftypes[af_idx]  #This is the array of airfoil objects that correspond to the airfoil index values

# Outputs before varying anything, ie "Normal"
Normalsections = Section.(r, chord, twist, airfoils)    #This is the reference with everything unchanged
Normalop = simple_op.(Vinf, Omega, r, rho) #The Normal operating points
Normalout = solve.(Ref(rotor), Normalsections, Normalop)  #outputs a struct of results for each section/radial location

#Vary the Variables- Rtip, chord, twist
```
Vary the radius- Note, r is a function of Rtip, so first take that out before redefining r
```
#Input Prep 
NumofRadiis = NumofVars #Number of radii to vary
NewRTip = range(5/100, 25/100, length = NumofRadiis) #This is the new Rtip values to vary over
rvar = r/Rtip #making rvar a function of NewRtip so Rtip can be varied 
chordrvar = chord/Rtip #making the chord a function of NewRtip so Rtip can be varied
EffRVar = zeros(NumofRadiis) #Initialize the efficiency array
CTRVar = zeros(NumofRadiis) #Initialize the CT array
CQRVar = zeros(NumofRadiis) #Initialize the CQ array

#Calculate and output
for i = 1:NumofRadiis
    if rvar[1].*NewRTip[i] < Rhub
        NewRHub = rvar[1].*NewRTip[i] #This is the new hub radius
        #println("Rhub was updated to $NewRHub for RVar")
    else
        NewRHub = Rhub
    end
    RVarRotor = Rotor(NewRHub, NewRTip[i], NumofBlades) #This is the new rotor object
    sections = Section.(rvar.*NewRTip[i], chordrvar.*NewRTip[i], twist, airfoils)
    op = simple_op.(Vinf, Omega, rvar.*NewRTip[i], rho) 
    out = solve.(Ref(RVarRotor), sections, op)
    T, Q = thrusttorque(RVarRotor, sections, out)   #calcs T & Q at each sec w/given conds, sums them for the whole rotor
    EffRVar[i], CTRVar[i], CQRVar[i] = nondim(T, Q, Vinf, Omega, rho, RVarRotor, "helicopter")
end
#Plotting
RVarCoefPlot = Plots.plot(NewRTip, CQRVar, title = "Varying Blade Radius", label = "Coefficient of Torque", legend = :topleft, margin = pad, box = :on,)
plot!(NewRTip, CQRVar.*2 .*pi, label = "Coefficient of Power", linecolor = :orange, ylabel = "Coef of Torque and Power",xlabel = "Radius (m)")
plot!(Plots.twinx(), NewRTip, CTRVar, label = "Coefficient of Thrust", legend = (.6, .7), linecolor = :green, ylabel = "Coef of Thrust")
RVarEffPlot = Plots.plot(NewRTip, EffRVar, label = "Figure of Merit", legend = :topright, margin = pad, title = "Varying Blade Radius", xlabel = "Radius (m)")
RVarPlot = Plots.plot(RVarCoefPlot, RVarEffPlot, layout = dims)

#Redo the above but vary RPM and plot as a surface
RPMVar = range(1000, 10000, length = NumofVars) #This is the RPM values to vary over
#? how do I deal with RHub as I can't run an if statement
function RSurfaceFunc(SRPM,SNewRTip;spitout) #! I think I need to do a nested for loop w/rpm and radius
    println(SRPM)
    println(SNewRTip)
    for i = 1:NumofRadiis
        if rvar[1]*SNewRTip[i] < Rhub
            NewRHub = rvar[1]*SNewRTip[i] #This is the new hub radius
            #println("Rhub was updated to $NewRHub for RVar")
        else
            NewRHub = Rhub
        end
        for k = 1:NumofVars
            SOmega = (pi/30)*SRPM[k]
            RVarRotor = Rotor(NewRHub, SNewRTip[i], NumofBlades) #This is the new rotor object
            sections = Section.(rvar.*SNewRTip[i], chordrvar.*SNewRTip[i], twist, airfoils)
            op = simple_op.(Vinf, SOmega, rvar.*SNewRTip[i], rho) 
            out = solve.(Ref(RVarRotor), sections, op)
            T, Q = thrusttorque(RVarRotor, sections, out)   #calcs T & Q at each sec w/given conds, sums them for the whole rotor
            EffRVar[k], CTRVar[k], CQRVar[k] = nondim(T, Q, Vinf, Omega, rho, RVarRotor, "helicopter")
            #println(EffRVar)
        end
    end
    if spitout == "Eff"
        println(EffRVar)
        return EffRVar
    elseif spitout == "CT"
        return CTRVar
    elseif spitout == "CQ"
        return CQRVar
    else
        println("Keyword argument error")
    end

end
println("did we get here") #!the above function is broken and shouldn't be trusted
RVarSurfacePlot = surface(RPMVar, NewRTip, RSurfaceFunc(RPMVar, NewRTip; spitout = "Eff"), title = "Varying RPM and Blade Radius", label = "Efficiency", legend = :topleft, margin = pad, box = :on)


```
Vary the chord- Because we are varying this independently of Rtip we dont need to divide it out,
    however, we will anyway as varying the chord as a percent of Rtip is more intuitive and robust
    than varying by fixed amounts.
```
#Input Prep
CVarchord = chord/Rtip #this makes it so we can vary chord as a percent of Rtip, rather than a fixed measurement
NumofChords = NumofVars #Number of chords to vary
VarChordPercent = range(-.09, .09, length = NumofChords) #The OG min % is .097 (so we can't subtract more than that) max % is .28
EffCVar = zeros(NumofChords) #Initialize the efficiency array
CTCVar = zeros(NumofChords) #Initialize the CT array
CQCVar = zeros(NumofChords) #Initialize the CQ array

#Calculate and output
for i = 1:NumofChords
    sections = Section.(r, ((CVarchord.+VarChordPercent[i]).*Rtip), twist, airfoils)
    op = simple_op.(Vinf, Omega, r, rho)
    out = solve.(Ref(rotor), sections, op)
    T, Q = thrusttorque(rotor, sections, out)   #calcs T & Q at each sec w/given conds, sums them for the whole rotor
    EffCVar[i], CTCVar[i], CQCVar[i] = nondim(T, Q, Vinf, Omega, rho, rotor, "helicopter")
    # calcs the coef of the blade under the given conditions at each advance ratio
end
#Plotting
CVarCoefPlot = Plots.plot(VarChordPercent, CTCVar, title = "Varying Blade Chord", label = "Coefficient of Thrust", legend = :topleft,xlabel = "Values used to vary Chord (in % of RTip)", margin = pad, ylabel = "Coef of Thrust")
plot!(Plots.twinx(), VarChordPercent, CQCVar, label = "Coefficient of Torque", legend = (.56, .44),linecolor = :orange, ylabel = "Coef of Torque")
CVarEffPlot = Plots.plot(VarChordPercent, EffCVar, label = "Figure of Merit", legend = :topright, margin = pad, title = "Varying Blade Chord",xlabel = "Values used to vary Chord (in % of RTip)")
CVarPlot = Plots.plot(CVarCoefPlot, CVarEffPlot, layout = dims)


```
Vary the Pitch- Because we are varying this independently of Rtip we dont need to divide it out
    Pitch is a seperate input argument in the op function, so the sections can be defined outside
    of the loop. This will "Pitch" the whole blade up or down, think like angle of attack 
``` 
#Input Prep
NumofPitches = NumofVars #Number of pitches to vary
VarPitch = range(-2, 25, length = NumofPitches) #Ranging from pitching down to pitching up
#! Pitch has to be limited to about -5 deg or Thrust goes neg and nondim throws an error
EffPVar = zeros(NumofPitches) #Initialize the efficiency array
CTPVar = zeros(NumofPitches) #Initialize the CT array
CQPVar = zeros(NumofPitches) #Initialize the CQ array
PVarSections = Section.(r, chord, twist, airfoils) #Defining the sections to use

#Calculate and output
for i = 1:NumofPitches
    pitch = VarPitch[i]*pi/180 #This is the pitch angle to use
    op = simple_op.(Vinf, Omega, r, rho; pitch)
    out = solve.(Ref(rotor), PVarSections, op)
    T, Q = thrusttorque(rotor, PVarSections, out)   #calcs T & Q at each sec w/given conds, sums them for the whole rotor
    EffPVar[i], CTPVar[i], CQPVar[i] = nondim(T, Q, Vinf, Omega, rho, rotor, "helicopter")
    # calcs the coef of the blade under the given conditions at each advance ratio
end
#Plotting
PVarCoefPlot = Plots.plot(VarPitch, CTPVar, title = "Varying Blade Pitch", label = "Coefficient of Thrust", legend = :topleft, xlabel = "Pitch added or subtracted from current value (deg)", margin = pad, ylabel = "Coef of Thrust")
plot!(Plots.twinx(), VarPitch, CQPVar, label = "Coefficient of Torque", legend = (.6,.4), linecolor = :orange, ylabel = "Coef of Torque")
PVarEffPlot = Plots.plot(VarPitch, EffPVar, label = "Figure of Merit", legend = :topright, title = "Varying Blade Pitch", margin = pad, xlabel = "Pitch added or subtracted from current value (deg)")
PVarPlot = Plots.plot(PVarCoefPlot, PVarEffPlot, layout = dims, xlabel = "Pitch added or subtracted from current value (deg)")


```
Vary the advance ratios- This is just like in the BEM Project, see BEM Project folder and Jupyter NB
```
#Input Prep
NumofAdvanceRatios = NumofVars #Number of advance ratios to vary
VarAdvanceRatio = range(0.1, .6, length = NumofAdvanceRatios) #The forward speed  goes from 1/10th of rpm to 60% rpm
#! J has to be limited to .6 as above that Thrust goes neg (drag) and nondim throws an error 
n = Omega/(2*pi)    # converts radians persecond to just rotations per second, same as rpm/60
D = 2 * Rtip    #Diameter of the prop is 2* the radius (note: Rtip is defined from center of rotation so includes Rhub)
EffJVar = zeros(NumofAdvanceRatios) #Initialize the efficiency array
CTJVar = zeros(NumofAdvanceRatios) #Initialize the CT array
CQJVar = zeros(NumofAdvanceRatios) #Initialize the CQ array

#Calculate and output
for i = 1:NumofAdvanceRatios
    local Vinf = VarAdvanceRatio[i] * D * n   #makes a local inflow veloc var at each advance ratio 
    local op = simple_op.(Vinf, Omega, r, rho)  #creates op pts at each blade section/location
    outputs = solve.(Ref(rotor), Normalsections, op) #uses all data from above plus local op conditions
    T, Q = thrusttorque(rotor, Normalsections, outputs)   #calcs T & Q at each sec w/given conds, sums them for the whole rotor
    EffJVar[i], CTJVar[i], CQJVar[i] = nondim(T, Q, Vinf, Omega, rho, rotor, "helicopter")
    # calcs the coef of the blade under the given conditions at each advance ratio
end
#Plotting
JVarCoefPlot = Plots.plot(VarAdvanceRatio, CTJVar, title = "Varying Advance Ratio", label = "Coefficient of Thrust", legend = (.001, .45), xaxis = "Advance Ratio (J)", margin = pad, ylabel = "Coef of Thrust")
plot!(Plots.twinx(), VarAdvanceRatio, CQJVar, label = "Coefficient of Torque", legend = :topright, linecolor = :orange, ylabel = "Coef of Torque")
JVarEffPlot = Plots.plot(VarAdvanceRatio, EffJVar, label = "Figure of Merit", legend = :topright, margin = pad, title = "Varying Advance Ratio", xaxis = "Advance Ratio (J)")
JVarPlot = Plots.plot(JVarCoefPlot, JVarEffPlot, layout = dims)


#? Do I want to vary the twist angle also? what would varying twist look like?

#Plot Palooza
#plot(RVarPlot, CVarPlot, PVarPlot, JVarPlot, layout = (2, 2))
#plot!(size = (1800,1000))

#plot(RVarCoefPlot,RVarEffPlot,CVarCoefPlot,CVarEffPlot,PVarCoefPlot,PVarEffPlot,JVarCoefPlot,JVarEffPlot,layout = (2,4))
#=
PlotArray = []
push!(PlotArray, RVarCoefPlot)
push!(PlotArray, RVarEffPlot)
push!(PlotArray, CVarCoefPlot)
push!(PlotArray, CVarEffPlot)
push!(PlotArray, PVarCoefPlot)
push!(PlotArray, PVarEffPlot)
push!(PlotArray, JVarCoefPlot)
push!(PlotArray, JVarEffPlot)
plot(PlotArray...)
#plot!(size = (3500,800))    #changes the plot size 
=#
display(RVarCoefPlot)
display(RVarEffPlot)
display(CVarCoefPlot)
display(CVarEffPlot)
display(PVarCoefPlot)
display(PVarEffPlot)
display(JVarCoefPlot)
display(JVarEffPlot)
display(RVarSurfacePlot)
