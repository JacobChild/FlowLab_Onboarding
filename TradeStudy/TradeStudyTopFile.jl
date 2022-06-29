#TradeStudyTopFile.jl
#Jacob Child
#June 28th, 2022
#PseudoCode: To practice reading files
# cd("C:/Users/child/Documents/Flow_Lab/Onboarding/TradeStudy")

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
aftypes = Array{AlphaAF}(undef, NumofSections) #Initialize the array of alphaAF objects

#CIncludes needed files for the functions
include("TSPolarPlotterFunction.jl") 
include("TSExtenderNSmootherFunction.jl")

#For loop to create the airfoil data and sections
for i in 1:NumofSections
    alphadegs, Cl, Cd, Cdp, Cm = TSPolarPlotterFunction(ContourX[i], ContourY[i], Re, M, Ncrit, alfa)
    alpharads = alphadegs * pi/180
    OldPlotArray, newalphadegs, ExtendedSmoothCl, ExtendedSmoothCd = TSExtenderNSmootherFunction(alpharads, Cl, Cd, Cdp, Cm, cr75)
    aftypes[i] = AlphaAF("CurrentAirfoilData.txt");   #this gives angle of attack, lift coefficient, drag coefficient
end

af_idx = zero(r)
#For loop to put which airfoil is used for each section into the af_idx array
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
sections = Section.(r, chord, twist, airfoils)