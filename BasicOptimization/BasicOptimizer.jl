#= 
BasicOptimizer.jl
Jacob Child
August 6, 2022

Pseudocode- Run all of the airfoil setup code so that all of that isn't rerun every time. Call an
constraint function that calculates torque and keeps it within bounds. Then Call the
objective function, that will take in an x in this format x = [x1, x2] and returns just Thrust.


    cd("C:/Users/child/Documents/Flow_Lab/Onboarding/BasicOptimization")
=#
using Evolutionary, Revise, Xfoil, DataFrames, FLOWMath, Dierckx, DelimitedFiles, Plots.Measures, Plots, CCBlade
include("OptimizationFunctionFile.jl")

#Plotting Settings
dims = (1,2) #for plot layouts after each parameter that was varied, ie "subplot" layouts
pad = 16mm #padding/margin for the plots

#File Parameters
NumofVars = 50 #length of the variable vectors, ie if =20 then 20 pitches will be tested
RPMVar = range(1000, 10000, length = NumofVars) #This is the RPM values to vary over
#Fluid Conditions 
rho = 1.225 #kg/m^3
mu = 1.78E-5 #kg/m/s
Vinf = 0  #inflow velocity, this means that the solve done below is for a singular advance ratio (hover in this case)
RPM = 5000
Omega = RPM*pi/30  #converts rpm to rad/s, derivation: rpm*pi*360deg/(60sec*180deg)-> rpm*pi/30
#Find where RPMVar is closest to RPM 
RPMLocation = FindFunc(RPMVar, RPM)[1] #finds the location of RPM in the RPMVar vector

#Import the high level csv file
CmdFileName = "DJI-II.csv"
cmdfile = readdlm("data\\rotors\\$CmdFileName", ',')
Rtip = cmdfile[2,2]
#Rhub = 0.004938276 #This is the first r value in the DJI-II prop, there is an interference, but that is fixed below
Rhub = cmdfile[3,2]
NumofBlades = cmdfile[4,2]
BladeFileName = cmdfile[5,2]

#Define the rotor object
MyRotor= Rotor(Rhub, Rtip, NumofBlades) #stores a "Rotor" object

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
    MyRotor= Rotor(Rhub, Rtip, NumofBlades) #redefines and stores the "Rotor" object
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


#Includes needed files for the functions
include("OptimizerPolarPlotterFunction.jl") 
include("OptimizerExtenderNSmootherFunction.jl")
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
    #OldPlotArray, newalphadegs, ExtendedSmoothCl, ExtendedSmoothCd = TSExtenderNSmootherFunction(alpharads, Cl, Cd, Cdp, Cm, cr75, Re, M)
    newalphadegs, ExtendedSmoothCl, ExtendedSmoothCd = TSExtenderNSmootherFunction(alpharads, Cl, Cd, Cdp, Cm, cr75, Re, M)
    aftypes[i] = AlphaAF("CurrentAirfoilData.txt");   #this gives angle of attack, lift coefficient, drag coefficient
end


#For loop to put which airfoil is used for each section into the af_idx array
af_idx = zero(r) #? Why does zero(r) seem to work in this case, but FMRVar = zero(NumofRadiis) doesn't work
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
Normalout = solve.(Ref(MyRotor), Normalsections, Normalop)  #outputs a struct of results for each section/radial location

#Calculation setup
CVarchord = chord/Rtip ##this makes it so we can vary chord as a percent of Rtip, rather than a fixed measurement


#Optimization setup
#Bounds 
lower = [-.075, -3] # in the format [chord, pitch]
upper = [.025, 8] # in the format [chord, pitch]
lc = [-Inf] #Lower constraint on Torque
uc = [.06] #Upper constraint on Torque, the old method gave .05998
constraints = WorstFitnessConstraints(lower, upper, lc, uc, ConstraintFunction)

#Run the Optimizer 
result = Evolutionary.optimize(ObjectiveFunction, GA(populationSize = 50)) #without constraints
#result = Evolutionary.optimize(ObjectiveFunction, constraints, GA(populationSize = 50))
