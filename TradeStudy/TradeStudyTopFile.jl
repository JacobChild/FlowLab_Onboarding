#TradeStudyTopFile.jl
#Jacob Child
#June 28th, 2022
#PseudoCode: To practice reading files
# cd("C:/Users/child/Documents/Flow_Lab/Onboarding/TradeStudy")

#Open and import csv file 
using CSV 

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
Locations = AirfoilFiles["r/R"]
NumofLocations = length(Locations)
ContourFiles = AirfoilFiles["Contour file"] #These are the files that contain the airfoil coordinates
AeroFiles = AirfoilFiles["Aero file"]   #These files have the xfoil output etc, but aren't extrapolated
#AeroFiles may not exist depending on the csv file and whether or not the data was previously calculated
#I will ignore the AeroFiles data in favor of more robust code that will generate my own

#Import the ContourFiles
ContourX = [Vector{Float64}(undef,(0)) for _ in 1:length(ContourFiles)] #Initialize the array of vectors
ContourY = [Vector{Float64}(undef,(0)) for _ in 1:length(ContourFiles)] #Initialize the array of vectors

#For loop that extracts the x and y coordinates from each of the ContourFiles
for i = 1:NumofLocations
    println(ContourFiles[i])
    GeomFileName = ContourFiles[i]
    GeomFiles = CSV.File("data\\airfoils\\$GeomFileName",delim = ",")
    ContourX[i] = GeomFiles["x/c"]
    ContourY[i] = GeomFiles["y/c"]
end

