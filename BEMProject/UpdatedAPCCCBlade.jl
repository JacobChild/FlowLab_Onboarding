#=
3DRCCBladeQuickStart.jl
Jacob Child 
June 1, 2022

Following the tutorial - https://flow.byu.edu/CCBlade.jl/stable/tutorial/
Goal: Follow the tutorial, but using an APC propeller 10x7 Thin Electric
The propeller is sold here: https://www.apcprop.com/product/10x7e/
And data is here: https://m-selig.ae.illinois.edu/props/volume-1/propDB-volume-1.html#APC
    Found under "Thin Electric" it is the "10 x 7  [1]" ctrl+f for that exact and you will find it under apc thin electric
10 x 7 means 10 in diameter and 7 inches pitch

PseudoCode- from UIUC do an analysis of a thin 10 x 7 propeller, 10 means diameter of 10inches,
    7 means 7in pitch, ie if screwed through a viscous material how far it would travel

cd("C:/Users/child/Documents/Flow_Lab/Onboarding/BEMProject")
=#

#Packages
using CCBlade, Plots

#Define the Rotor
Rtip = 10/2 * .0254     #The tip radius is the diameter/2 and the converts in to meters
Rhub = .8/2 * .0254     #The hub radiaus is diameter/2 and can be found on the website
#Note: if not given on a website, look at the geometry file and the first point is r/R, or saying the current point is 15% of the tip radius 
# ...  This means that the hub has to be less than that amount as the first point is on the blade itself and not the hub 
B = 2   #Number of blades
RPM = 5000.0  #The quick start guide uses 5400, but experimental data is only in 1000s on the website 

rotor = Rotor(Rhub, Rtip, B) #calls Rotor from the CCBlade package and stores a "Rotor" object?
#fieldnames(typeof(rotor))

#Define the Prop geometry
#Taken from the UIUC database for this propDB
#The first 2 columns are given in percent of Rtip (or Blade Tip Radius)

propgeom = [
0.15   0.138   37.86
0.20   0.154   45.82
0.25   0.175   44.19
0.30   0.190   38.35
0.35   0.198   33.64
0.40   0.202   29.90
0.45   0.200   27.02
0.50   0.195   24.67
0.55   0.186   22.62
0.60   0.174   20.88
0.65   0.161   19.36
0.70   0.145   17.98
0.75   0.129   16.74
0.80   0.112   15.79
0.85   0.096   14.64
0.90   0.081   13.86
0.95   0.061   12.72
1.00   0.040   11.53
]

r = propgeom[:, 1] * Rtip   #This gives the radial locations r in 5% steps along the blade
chord = propgeom[:, 2] * Rtip   #This gives the Chord length at each radial step/location 
theta = propgeom[:, 3] * pi/180     #the twist angle, converted from degrees to radians (what CCBlade wants)


#Airfoil data
#Airfoil data needs to be imported, ie the lift and drag coef of the airfoil at each location
#I haven't found the exact data needed, so we are assuming Naca 4412, as their website said it is close to that
af = AlphaAF("EpplerE63Data.txt")   #this gives angle of attack, lift coefficient, drag coefficient
#... the top row has header info, then reynolds number, then mach number, on a dif row each

#Define sections to analyze
sections = Section.(r, chord, theta, Ref(af))   #the "." broadcasts to define all sections at once, and as the airfoil 
#... doesn't change af is wrapped in Ref to make it the same for each 
#fieldnames(typeof(sections[1]))

#Define the operating points- or what the flow field/conditions are at each section 
#Define the input airflow and broadcast to solve the velocity at each point 
Vinf = 5.0  #inflow velocity
Omega = RPM*pi/30  #converts rpm to rad/s, derivation: rpm*pi*360deg/(60sec*180deg)-> rpm*pi/30
rho = 1.225     #density of the air in metric

op = simple_op.(Vinf, Omega, r, rho)    #calls simple_op to define "operating points" at each section
#fieldnames(typeof(op[1])) shows that each op pt contains the velocity components, density, pitch   
#mu (dynamic viscosity?), and the speed of sound (not sure why that is needed)


#Solve with the Blade Element Method 
out = solve.(Ref(rotor), sections, op)  #outputs a struct of results for each section/radial location?
# struct outputs can be found [here](https://flow.byu.edu/CCBlade.jl/stable/reference/#Output-Struct)

#Plots and Data analysis

#Plot Forces
#= figure()
plot(r/Rtip, out.Np)    #plots location (as percent of radius) and the normal force per unit length 
plot(r/Rtip, out.Tp) #plots "" and the tangential force per unit length 
xlabel("r/Rtip")
ylabel("distributed loads (N/m)")
legend(["flapwise (normal force", "lead-lag (tangential force)"])
 =#


#Plot induced velocities, can be used for future interaction calculations
#= figure()
plot(r/Rtip, out.u/Vinf) #axial induced velocity/Inflow velocity, so it is a percent, ie normalized
plot(r/Rtip, out.v/Vinf) #tangential induced velocity (swirl)/""...
xlabel("r/Rtip")
ylabel("Induced velocity at rotor disk (normalized)")
legend(["axial velocity", "swirl velocity"]) =#


#Prep for Coefficient analysis (thrust, power efficiency etc) at different advance ratios
#Advance ratio = ratio of forward speed/ (rotational speed times prop diameter)

#Inputs
nJ = 30     #number of advance ratios to evaluate
J = range(.00001, .92, length = nJ) #creates a range of advance ratios from .001 to .92 (exp ending), what are normal vals?
Omega = RPM*pi/30  #rpms to rad/sec 
n = Omega/(2*pi)    # converts radians persecond to just rotations per second, same as rpm/60
D = 2 * Rtip    #Diameter of the prop is 2* the radius, do I ignore the hub?

eff = zeros(nJ)     #creates arrays for efficiency, coef of thrust, and coef of torque 
CT = zeros(nJ)
CQ = zeros(nJ)      #coef of torque is requried torque over theoretical required torque
T = zeros(nJ)
Q = zeros(nJ) #to define them outside of the for loop

#Coef Solver

for i = 1:nJ
    local Vinf = J[i] * D * n   #makes a local inflow veloc var at each advance ratio 
    local op = simple_op.(Vinf, Omega, r, rho)  #creates op pts at each blade section/location

    outputs = solve.(Ref(rotor), sections, op) #uses all data from above plus local op conditions
    T[i], Q[i] = thrusttorque(rotor, sections, outputs)   #calcs T & Q at each sec w/given conds, sums them for the whole rotor
    eff[i], CT[i], CQ[i] = nondim(T[i], Q[i], Vinf, Omega, rho, rotor, "propeller")
    # calcs the coef of the blade under the given conditions at each advance ratio
end

#Prep to compare to experimental data from https://www.apcprop.com/files/PER3_10x7E.dat
#Note: this experimental data is specific to this prop at 5000rpm
#Headers: V (mph) | J (Adv Ratio) | Pe (efficiency) | Ct | Cp | PWR (Hp) | Torque (In-Lbf) | Thrust (Lbf)

exp = [
    0.0        0.00      0.0000      0.1161      0.0549       0.055       0.695       0.925              
    1.5        0.03      0.0657      0.1159      0.0559       0.056       0.708       0.923              
    3.0        0.06      0.1286      0.1157      0.0569       0.057       0.722       0.921              
    4.5        0.09      0.1888      0.1153      0.0580       0.058       0.735       0.918              
    6.0        0.13      0.2462      0.1149      0.0591       0.059       0.749       0.915              
    7.5        0.16      0.3008      0.1143      0.0601       0.060       0.762       0.910              
    9.0        0.19      0.3527      0.1133      0.0610       0.061       0.774       0.902              
   10.5        0.22      0.4018      0.1119      0.0617       0.062       0.782       0.891              
   12.0        0.25      0.4481      0.1100      0.0622       0.063       0.788       0.876              
   13.5        0.28      0.4915      0.1077      0.0624       0.063       0.791       0.858              
   15.0        0.32      0.5318      0.1050      0.0625       0.063       0.792       0.836              
   16.5        0.35      0.5691      0.1018      0.0623       0.063       0.790       0.811              
   18.0        0.38      0.6033      0.0983      0.0619       0.062       0.785       0.783              
   19.5        0.41      0.6343      0.0945      0.0613       0.062       0.777       0.752              
   21.0        0.44      0.6620      0.0902      0.0604       0.061       0.765       0.718              
   22.5        0.47      0.6864      0.0855      0.0591       0.059       0.750       0.681              
   24.0        0.51      0.7076      0.0803      0.0575       0.058       0.729       0.640              
   25.5        0.54      0.7256      0.0749      0.0556       0.056       0.704       0.597              
   27.0        0.57      0.7406      0.0694      0.0534       0.054       0.677       0.553              
   28.5        0.60      0.7526      0.0637      0.0509       0.051       0.646       0.508              
   30.0        0.63      0.7624      0.0580      0.0482       0.048       0.610       0.462              
   31.5        0.66      0.7701      0.0521      0.0450       0.045       0.570       0.415              
   33.0        0.70      0.7755      0.0460      0.0413       0.042       0.524       0.366              
   34.5        0.73      0.7781      0.0398      0.0372       0.037       0.472       0.317              
   36.0        0.76      0.7773      0.0334      0.0327       0.033       0.414       0.266              
   37.5        0.79      0.7694      0.0269      0.0277       0.028       0.351       0.214              
   39.0        0.82      0.7468      0.0203      0.0224       0.023       0.284       0.162              
   40.5        0.85      0.6956      0.0136      0.0167       0.017       0.212       0.108              
   42.0        0.89      0.5648      0.0068      0.0107       0.011       0.136       0.054              
   43.5        0.92      0.0024      0.0000      0.0047       0.005       0.059       0.000 
] 
#Above is the experimental data, I assume at each radial location 

#Experimental values
JExp = exp[:,2]     #Advanced ratio
CtExp = exp[:, 4] #Coef of Thrust
CpExp = exp[:, 5] #Coef of power 
EtaExp = exp[:, 3] #efficiency


#= figure()
plot(J, CT)
plot(J, CQ*2*pi) #converts coef of torque to power
plot(JExp, CtExp, "ko")
plot(JExp, CpExp, "ko")
xlabel(L"J")
legend([L"C_T", L"C_P", "experimental"])

figure()
plot(J, eff)
plot(JExp, EtaExp, "ko")
xlabel(L"J")
ylabel(L"\eta")
legend(["CCBlade", "experimental"]) =#

default(show=true)

CTPlot = plot(J, CT, title = "Coefficient of Thrust", label = "Theoretical Ct")
ExpCTPlot = plot!(JExp, CtExp, m = 4, xlabel = "Advance Ratio (J)", label = "Experimental");

#Error Calculations
RelError = abs.(CT - CtExp) ./ CtExp * 100
CTErrorPlot = plot(J,RelError, title = "Relative Error", xlabel = "Advance Ratio (J)", ylabel = "Percent Relative Error",
    label = "Error %")

#plot(ExpCTPlot,CTErrorPlot)

CPPlot = plot(J, CQ*2*pi, title = "Coefficient of Power", label = "Theoretical Cp") #CQ*2*pi converts coef of torque to coef of power
ExpCPPlot  = plot!(JExp, CpExp, m = 4, xlabel = "Advance Ratio (J)", label = "Experimental");

#Error Calculations
RelError = (CQ*2*pi - CpExp) ./ CpExp * 100
CPErrorPlot = plot(J,RelError, title = "Relative Error", xlabel = "Advance Ratio (J)", 
ylabel = "Percent Relative Error", label = "Error %")

#plot(ExpCPPlot, CPErrorPlot)

EffPlot = plot(J[1:29], eff[1:29], title = "Effeciency", label = "Theoretical Eta")
ExpEffPlot = plot!(JExp[1:29], EtaExp[1:29], m = 4, xlabel = "Advance Ratio (J)", ylabel = "Effeciency (eta)",
    label = "Experimental", legend = :topleft);

#Error Calculations
RelError = (eff[1:29] - EtaExp[1:29]) ./ EtaExp[1:29] * 100
EffErrorPlot = plot(J[1:29],RelError, title = "Relative Error", xlabel = "Advance Ratio (J)", 
ylabel = "Percent Relative Error", label = "Error %")

#plot(ExpEffPlot,EffErrorPlot)

#Torque and Thrust plots
TPlot = plot(J, T, label = "Thrust (T)")
QPlot = plot(J, Q, label = "Torque(Q)")
plot(TPlot, QPlot, ExpCTPlot, ExpCPPlot, ExpEffPlot, legend = false)