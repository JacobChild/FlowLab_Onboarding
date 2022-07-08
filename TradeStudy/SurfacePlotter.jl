#=
SurfacePlotter.jl
Jacob Child
July 1st, 2022

Purpose: To learn how to plot surfaces so I can better understand my data and implement surface 
plotting in TradeStudyTopFile.

Pseudocode: Copy and import all the needed data for varying the pitch, from the TradeStudyTopFile.jl
and from the REPL, once that is hard coded in, copy the old for loops and turn them into a function
that can be evaluated at a range of RPMs and range of pitches. The output can be Coef of Thrust or 
efficiency.

cd("C:/Users/child/Documents/Flow_Lab/Onboarding/TradeStudy")

=#



#Function Fun 
function CoefThrustFunc(SRPM,SVarPitch)
    #println(SRPM)
    SOmega = (pi/30).*SRPM  #converts rpm to rad/s, derivation: rpm*pi*360deg/(60sec*180deg)-> rpm*pi/30
    pitch = SVarPitch*pi/180 #This is the pitch angle to use
    Sop = simple_op.(Vinf, SOmega, Sr, rho; pitch)
    Sout = solve.(Ref(Srotor), SPVarSections, Sop)
    ST, SQ = thrusttorque(Srotor, SPVarSections, Sout)   #calcs T & Q at each sec w/given conds, sums them for the whole rotor
    #SEffPVar, SCTPVar, SCQPVar = nondim(ST, SQ, Vinf, SOmega, rho, Srotor, "propeller")
    SEffPVar, SCTPVar, SCQPVar = nondim(ST, SQ, Vinf, SOmega, rho, Srotor, "helicopter")
    #println( ST / (rho * pi * Rp^2 * (SOmega*Rp)^2))
    #println(SCTPVar)
    return SCTPVar
end

function CoefPowerFunc(SRPM,SVarPitch)
SOmega = (pi/30).*SRPM  #converts rpm to rad/s, derivation: rpm*pi*360deg/(60sec*180deg)-> rpm*pi/30
#local Omega = SOmega
#println(SOmega)
pitch = SVarPitch*pi/180 #This is the pitch angle to use
Sop = simple_op.(Vinf, SOmega, Sr, rho; pitch)
Sout = solve.(Ref(Srotor), SPVarSections, Sop)
ST, SQ = thrusttorque(Srotor, SPVarSections, Sout)   #calcs T & Q at each sec w/given conds, sums them for the whole rotor
#SEffPVar, SCTPVar, SCQPVar = nondim(T, Q, Vinf, SOmega, rho, Srotor, "propeller")
SEffPVar, SCTPVar, SCQPVar = nondim(ST, SQ, Vinf, SOmega, rho, Srotor, "helicopter")

return SCQPVar
end

function EffFunc(SRPM,SVarPitch)
#println(Vinf)
SOmega = (pi/30).*SRPM  #converts rpm to rad/s, derivation: rpm*pi*360deg/(60sec*180deg)-> rpm*pi/30
#local Omega = SOmega
pitch = SVarPitch*pi/180 #This is the pitch angle to use
Sop = simple_op.(Vinf, SOmega, Sr, rho; pitch)
Sout = solve.(Ref(Srotor), SPVarSections, Sop)
ST, SQ = thrusttorque(Srotor, SPVarSections, Sout)   #calcs T & Q at each sec w/given conds, sums them for the whole rotor
#SEffPVar, SCTPVar, SCQPVar = nondim(ST, SQ, Vinf, SOmega, rho, Srotor, "propeller")
SEffPVar, SCTPVar, SCQPVar = nondim(ST, SQ, Vinf, SOmega, rho, Srotor, "helicopter")
# calcs the coef of the blade under the given conditions at each advance ratio

return SEffPVar
end

function ThrustFunc(SRPM,SVarPitch)
    #println(Vinf)
    SOmega = (pi/30).*SRPM  #converts rpm to rad/s, derivation: rpm*pi*360deg/(60sec*180deg)-> rpm*pi/30
    #local Omega = SOmega
    pitch = SVarPitch*pi/180 #This is the pitch angle to use
    Sop = simple_op.(Vinf, SOmega, Sr, rho; pitch)
    Sout = solve.(Ref(Srotor), SPVarSections, Sop)
    ST2, SQ = thrusttorque(Srotor, SPVarSections, Sout)   #calcs T & Q at each sec w/given conds, sums them for the whole rotor
    #println(ST2)
    return ST2
end


#Needed data:
Vinf = 0
Sr = [0.004938276, 0.008230452, 0.011440332, 0.0141564, 0.0165432, 0.01868316, 0.02115228, 0.02403288, 0.02683128, 0.030041159999999997, 0.0330864, 0.03604944, 0.03893004, 0.04197528, 0.04543212, 0.04987656, 0.054650159999999996, 0.059670839999999996, 0.06395064, 0.06921816, 0.07333332, 0.07827156, 0.0826338, 0.09201648, 0.10131684, 0.11061732, 0.12]
rho = 1.225
SNumofPitches = 50 #Number of pitches to vary
#SVarPitch = fill(10,SNumofPitches)
SVarPitch = range(0, 25, length = SNumofPitches) #Ranging from pitching down to pitching up
SRPM = range(3000, 10000, length = SNumofPitches) #Ranging from 1000 to 10000 RPMs
#! Omega is a Function of RPM Omega = RPM*pi/30  
Schord = [0.014521319999999999, 0.016580519999999998, 0.01819272, 0.01990584, 0.021340199999999997, 0.022876439999999998, 0.02493072, 0.027620759999999998, 0.030347639999999995, 0.032067479999999995, 0.03361188, 0.033866759999999996, 0.03303648, 0.03160032, 0.030318239999999996, 0.028672919999999998, 0.026818679999999998, 0.02521956, 0.024161759999999997, 0.02311272, 0.02220264, 0.02101788, 0.0202266, 0.01801932, 0.01602948, 0.013877039999999998, 0.011740332]
Stwist = [0.28722359901295086, 0.30543261909900765, 0.323186108250294, 0.3382081571222092, 0.3514098275842943, 0.3632449052420678, 0.3769003613096714, 0.39283347205112773, 0.4083110518578134, 0.42606454100909974, 0.44290696829084497, 0.45851893844993435, 0.4314366644467382, 0.4068467206153902, 0.3829287285460599, 0.35704898639748794, 0.3339408271010831, 0.31362693993712104, 0.29882829320946114, 0.28312731125852014, 0.2724301882730469, 0.2610785668180758, 0.2521791329621567, 0.2358969563703016, 0.22273368315176034, 0.21178523275399994, 0.20245819323134223]
Sairfoils = airfoils; #this was too large to see and copy in the terminal, so it was just assigned as it is in the REPL memory
#! The above won't work if you haven't run the TradeStudyTopFile.jl yet
Srotor = Rotor{Float64, Int64, Bool, Nothing, Nothing, Nothing, PrandtlTipHub}(0.004938276, 0.12, 2, 0.0, false, nothing, nothing, nothing, PrandtlTipHub())
SPVarSections = Section.(Sr, Schord, Stwist, Sairfoils) #Defining the sections to use

#Plotting
using PyCall #allows figures you can interact with, ie rotate
pygui(true) #makes it so I can rotate the plots
using PyPlot
CoefofThrustPlot = surface(SRPM, SVarPitch, CoefThrustFunc, xlabel = "RPM",ylabel = "Pitch (deg)", zlabel = "Coef of Thrust")
CoefofEffPlot = surface(SRPM, SVarPitch, EffFunc, xlabel = "RPM",ylabel = "Pitch (deg)", zlabel = "Efficiency")
CoefofPowerPlot = surface(SRPM, SVarPitch, CoefPowerFunc, xlabel = "RPM",ylabel = "Pitch (deg)", zlabel = "Power")
ThrustSurfacePlot = surface(SRPM, SVarPitch, ThrustFunc, xlabel = "RPM",ylabel = "Pitch (deg)", zlabel = "Thrust")

display(CoefofThrustPlot)
display(CoefofEffPlot)
display(CoefofPowerPlot)
display(ThrustSurfacePlot)

SEffPVar = zeros(SNumofPitches)
SCTPVar = zeros(SNumofPitches)
SCQPVar = zeros(SNumofPitches)
ST1 = zeros(SNumofPitches)
SQ1 = zeros(SNumofPitches)

for i = 1:SNumofPitches

    SOmega = (pi/30).*SRPM[i]  #converts rpm to rad/s, derivation: rpm*pi*360deg/(60sec*180deg)-> rpm*pi/30
    pitch = SVarPitch[i]*pi/180 #This is the pitch angle to use
    Sop = simple_op.(Vinf, SOmega, Sr, rho; pitch)
    Sout = solve.(Ref(Srotor), SPVarSections, Sop)
    ST1[i], SQ1[i] = thrusttorque(Srotor, SPVarSections, Sout)   #calcs T & Q at each sec w/given conds, sums them for the whole rotor
    #SEffPVar, SCTPVar, SCQPVar = nondim(ST, SQ, Vinf, SOmega, rho, Srotor, "propeller")
    SEffPVar[i], SCTPVar[i], SCQPVar[i] = nondim(ST1[i], SQ1[i], Vinf, SOmega, rho, Srotor, "helicopter")
    #println( ST / (rho * pi * Rp^2 * (SOmega*Rp)^2))
    #println(SCTPVar)
end
StraightThrustPlot = plot(SRPM,ST1)