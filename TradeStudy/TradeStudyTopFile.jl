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
using CSV, Xfoil, CCBlade, Printf, DataFrames, Plots, FLOWMath, Dierckx, DelimitedFiles
#include("TSExtenderNSmootherFunction.jl")

#Open and import csv file
#Import the high level csv file
CmdFileName = "DJI-II.csv"
cmdfile = CSV.File("data\\rotors\\$CmdFileName",delim = ",")
Rtip = parse(Float64, cmdfile.file[1])
Rhub = parse(Float64, cmdfile.file[2])
NumofBlades = parse(Float64, cmdfile.file[3])
BladeFileName = cmdfile.file[4]

#Define the rotor object
rotor = Rotor(Rhub, Rtip, NumofBlades) #stores a "Rotor" object

#Import the blade directory file
BladeFile = CSV.File("data\\rotors\\$BladeFileName",delim = ",")
ChordFile = BladeFile.file[1]
PitchFile = BladeFile.file[2]
SweepFile = BladeFile.file[3]
HeightFile = BladeFile.file[4]
AirfoilFiles = BladeFile.file[5]
SplineOrder = parse(Int64, BladeFile.file[6])
SplineSmoothing = parse(Float64, BladeFile.file[7])

#Import the ChordFile 
ChordDist = CSV.File("data\\rotors\\$ChordFile",delim = ",")
r = ChordDist["r/R"] * Rtip
chord = ChordDist["c/R"] * Rtip

#Import the PitchFile
PitchDist = CSV.File("data\\rotors\\$PitchFile",delim = ",")
twist = PitchDist["twist (deg)"] * pi/180 #import twist and convert to radians

#Import the AirfoilFiles 
AirfoilFiles = CSV.File("data\\rotors\\$AirfoilFiles",delim = ",")
Locations = AirfoilFiles["r/R"]*Rtip
NumofSections = length(Locations)
ContourFiles = AirfoilFiles["Contour file"] #These are the files that contain the airfoil coordinates
AeroFiles = AirfoilFiles["Aero file"]   #These files have the xfoil output etc, but aren't extrapolated
#AeroFiles may not exist depending on the csv file and whether or not the data was previously calculated
#I will ignore the AeroFiles data in favor of more robust code that will generate my own

#Import the ContourFiles
ContourX = [Vector{Float64}(undef,(0)) for _ in 1:length(ContourFiles)] #Initialize the array of vectors
ContourY = [Vector{Float64}(undef,(0)) for _ in 1:length(ContourFiles)] #Initialize the array of vectors

#For loop that extracts the x and y coordinates from each of the ContourFiles
for i = 1:NumofSections
    GeomFileName = ContourFiles[i]
    GeomFiles = CSV.File("data\\airfoils\\$GeomFileName",delim = ",")
    ContourX[i] = GeomFiles["x/c"]
    ContourY[i] = GeomFiles["y/c"]
end

#Define the fluid conditions and other values for functions
#cr75
indexval = round(Int64, .75*length(chord))
cr75 = chord[indexval]/Rtip

#Fluid Conditions 
rho = 1.225 #kg/m^3
Re = 50000 #Reynolds number
M = 0.0     #Mach number
Ncrit = 9   #the "smoothness of the air"
alfa = -1.0:.5:10.0 #This is the angle of attack in degrees, any benefit to a larger range?
Vinf = 5.0  #inflow velocity, this means that the solve done below is for a singular advance ratio
RPM = 5000
Omega = RPM*pi/30  #converts rpm to rad/s, derivation: rpm*pi*360deg/(60sec*180deg)-> rpm*pi/30


#CIncludes needed files for the functions
include("TSPolarPlotterFunction.jl") 
include("TSExtenderNSmootherFunction.jl")

#For loop to create the airfoil data and sections
aftypes = Array{AlphaAF}(undef, NumofSections) #Initialize the array of alphaAF objects

for i in 1:NumofSections
    alphadegs, Cl, Cd, Cdp, Cm = TSPolarPlotterFunction(ContourX[i], ContourY[i], Re, M, Ncrit, alfa)
    alpharads = alphadegs * pi/180
    OldPlotArray, newalphadegs, ExtendedSmoothCl, ExtendedSmoothCd = TSExtenderNSmootherFunction(alpharads, Cl, Cd, Cdp, Cm, cr75)
    aftypes[i] = AlphaAF("CurrentAirfoilData.txt");   #this gives angle of attack, lift coefficient, drag coefficient
end


#For loop to put which airfoil is used for each section into the af_idx array
af_idx = zero(r)

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



#Vary the Variables- r, chord, twist 
sections = Section.(r, chord, twist, airfoils)    #This is the reference with everything unchanged
op = simple_op.(Vinf, Omega, r, rho) #The Normal operating points
#TODO Fix this from throwing an error, it does so in the next one, so I assume it won't work in all of them
#Normalout = solve.(Ref(rotor), sections, op)  #outputs a struct of results for each section/radial location

```
Vary the radius- Note, r is a function of Rtip, so first take that out before redefining r
```
#Input Prep 
NumofRadiis = 20 #Number of radii to vary
NewRTip = range(5, 25, length = NumofRadiis) #This is the new Rtip values to vary over
rvar = r/Rtip #making rvar a function of NewRtip so Rtip can be varied 
chordrvar = chord/Rtip #making the chord a function of NewRtip so Rtip can be varied
EffRVar = zero(NumofRadiis) #Initialize the efficiency array
CTRVar = zero(NumofRadiis) #Initialize the CT array
CQRVar = zero(NumofRadiis) #Initialize the CQ array

#Calculate and output
for i = 1:NumofRadiis
    sections = Section.(rvar*NewRTip[i], chordrvar*NewRTip[i], twist, airfoils)
    op = simple_op.(Vinf, Omega, rvar*NewRTip[i], rho)
    Normalout = solve.(Ref(rotor), sections, op)
    T, Q = thrusttorque(rotor, sections, Normalout)   #calcs T & Q at each sec w/given conds, sums them for the whole rotor
    EffRVar[i], CTRVar[i], CQRVar[i] = nondim(T, Q, Vinf, Omega, rho, rotor, "propeller")
    # calcs the coef of the blade under the given conditions at each advance ratio
end
RVarCoefPlot = plot(NewRtip, CTRVar, title = "Varying Blade Radius", label = "Coefficient of Thrust", legend = :topleft)
plot!(twinx(), NewRtip, CQRVar, label = "Coefficient of Torque", legend = :topright)
RVarEffPlot = plot(NewRtip, EffRVar, label = "Efficiency", legend = :topright)
RVarPlot = plot(RVarCoefPlot, RVarEffPlot, layout = (2,1))


```
Vary the chord- Because we are varying this independently of Rtip we don't need to divide it out,
    however, we will anyway as varying the chord as a percent of Rtip is more intuitive and robust
    than varying by fixed amounts.
```
#Input Prep
CVarchord = chord/Rtip #this makes it so we can vary chord as a percent of Rtip, rather than a fixed measurement
NumofChords = 20 #Number of chords to vary
VarChordPercent = range(-.09, .09, length = NumofChords) #The OG min % is .097 (so we can't subtract more than that) max % is .28
EffCVar = zero(NumofChords) #Initialize the efficiency array
CTCVar = zero(NumofChords) #Initialize the CT array
CQCVar = zero(NumofChords) #Initialize the CQ array

#Calculate and output
for i = 1:NumofChords
    sections = Section.(r, ((CVarchord.+VarChordPercent[i]).*Rtip), twist, airfoils)
    op = simple_op.(Vinf, Omega, r, rho)
    out = solve.(Ref(rotor), sections, op)
    T, Q = thrusttorque(rotor, sections, out)   #calcs T & Q at each sec w/given conds, sums them for the whole rotor
    EffCVar[i], CTCVar[i], CQCVar[i] = nondim(T, Q, Vinf, Omega, rho, rotor, "propeller")
    # calcs the coef of the blade under the given conditions at each advance ratio
end
CVarCoefPlot = plot(VarChordPercent, CTCVar, title = "Varying Blade Chord", label = "Coefficient of Thrust", legend = :topleft)
plot!(twinx(), VarChordPercent, CQCVar, label = "Coefficient of Torque", legend = :topright)
CVarEffPlot = plot(VarChordPercent, EffCVar, label = "Efficiency", legend = :topright)
CVarPlot = plot(CVarCoefPlot, CVarEffPlot, layout = (2,1))


```
Vary the Pitch- Because we are varying this independently of Rtip we don't need to divide it out
    Pitch is a seperate input argument in the op function, so the sections can be defined outside
    of the loop. This will "Pitch" the whole blade up or down, think like angle of attack 
``` 
#Input Prep
NumofPitches = 20 #Number of pitches to vary
VarPitch = range(-45, 45, Length = NumofPitches) #Ranging from pitching down to pitching up
EffPVar = zero(NumofPitches) #Initialize the efficiency array
CTPVar = zero(NumofPitches) #Initialize the CT array
CQPVar = zero(NumofPitches) #Initialize the CQ array
PVarSections = Section.(r, chord, twist, airfoils) #Defining the sections to use

#Calculate and output
for i = 1:NumofPitches
    op = simple_op.(Vinf, Omega, r, rho; VarPitch[i])
    out = solve.(Ref(rotor), PVarSections, op)
    T, Q = thrusttorque(rotor, PVarSections, out)   #calcs T & Q at each sec w/given conds, sums them for the whole rotor
    EffPVar[i], CTPVar[i], CQPVar[i] = nondim(T, Q, Vinf, Omega, rho, rotor, "propeller")
    # calcs the coef of the blade under the given conditions at each advance ratio
end
PVarCoefPlot = plot(VarPitch, CTPVar, title = "Varying Blade Pitch", label = "Coefficient of Thrust", legend = :topleft)
plot!(twinx(), VarPitch, CQPVar, label = "Coefficient of Torque", legend = :topright)
PVarEffPlot = plot(VarPitch, EffPVar, label = "Efficiency", legend = :topright)
PVarPlot = plot(PVarCoefPlot, PVarEffPlot, layout = (2,1))

```
Vary the advance ratios- This is just like in the BEM Project (see BEM Project folder and Jupyter NB)
```
#Input Prep
NumofAdvanceRatios = 20 #Number of advance ratios to vary
VarAdvanceRatio = range(0.1, .9, Length = NumofAdvanceRatios) #The forward speed  goes from 1/10th of rpm to 90% rpm 
n = Omega/(2*pi)    # converts radians persecond to just rotations per second, same as rpm/60
D = 2 * Rtip    #Diameter of the prop is 2* the radius (note: Rtip is defined from center of rotation so includes Rhub)
EffJVar = zero(NumofAdvanceRatios) #Initialize the efficiency array
CTJVar = zero(NumofAdvanceRatios) #Initialize the CT array
CQJVar = zero(NumofAdvanceRatios) #Initialize the CQ array

#Calculate and output
for i = 1:NumofAdvanceRatios
    local Vinf = VarAdvanceRatio[i] * D * n   #makes a local inflow veloc var at each advance ratio 
    local op = simple_op.(Vinf, Omega, r, rho)  #creates op pts at each blade section/location
    
    outputs = solve.(Ref(rotor), sections, op) #uses all data from above plus local op conditions
    T, Q = thrusttorque(rotor, sections, outputs)   #calcs T & Q at each sec w/given conds, sums them for the whole rotor
    EffJVar[i], CTJVar[i], CQJVar[i] = nondim(T, Q, Vinf, Omega, rho, rotor, "propeller")
    # calcs the coef of the blade under the given conditions at each advance ratio
end
JVarCoefPlot = plot(VarAdvanceRatio, CTJVar, title = "Varying Advance Ratio", label = "Coefficient of Thrust", legend = :topleft)
plot!(twinx(), VarAdvanceRatio, CQJVar, label = "Coefficient of Torque", legend = :topright)
JVarEffPlot = plot(VarAdvanceRatio, EffJVar, label = "Efficiency", legend = :topright)
JVarPlot = plot(JVarCoefPlot, JVarEffPlot, layout = (2,1))


#? Do I want to vary the twist angle also? what would varying twist look like?

#Plot
plot(RVarPlot, CVarPlot, PVarPlot, JVarPlot, layout = (4,1))