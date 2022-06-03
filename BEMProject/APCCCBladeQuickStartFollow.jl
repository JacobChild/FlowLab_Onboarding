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

PseudoCode- from UIUC do an analysis of a thin 10 x 5 propeller, 10 means diameter of 10inches

cd("C:/Users/child/Documents/Flow_Lab/Onboarding/BEMProject")
=#

#Packages
using CCBlade, PyPlot

#Define the Rotor
Rtip = 10/2 * .0254     #The tip radius is the diameter/2 and the converts in to meters
Rhub = .8/2 * .0254     #The hub radiaus is diameter/2 and can be found on the website
#Note: if not given on a website, look at the geometry file and the first point is r/R, or saying the current point is 15% of the tip radius 
# ...  This means that the hub has to be less than that amount as the first point is on the blade itself and not the hub 
B = 2   #Number of blades

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
af = AlphaAF("naca4412.dat")   #this gives angle of attack, lift coefficient, drag coefficient
#... the top row has header info, then reynolds number, then mach number, on a dif row each

#Define sections to analyze
sections = Section.(r, chord, theta, Ref(af))   #the "." broadcasts to define all sections at once, and as the airfoil 
#... doesn't change af is wrapped in Ref to make it the same for each 
#fieldnames(typeof(sections[1]))

#Define the operating points- or what the flow field/conditions are at each section 
#Define the input airflow and broadcast to solve the velocity at each point 
Vinf = 5.0  #inflow velocity
Omega = 5400*pi/30  #converts rpm to rad/s, derivation: rpm*pi*360deg/(60sec*180deg)-> rpm*pi/30
rho = 1.225     #density of the air in metric

op = simple_op.(Vinf, Omega, r, rho)    #calls simple_op to define "operating points" at each section
#fieldnames(typeof(op[1])) shows that each op pt contains the velocity components, density, pitch   
#mu (dynamic viscosity?), and the speed of sound (not sure why that is needed)


#Solve with the Blade Element Method 
out = solve.(Ref(rotor), sections, op)  #outputs a struct of results for each section/radial location?
# struct outputs can be found [here](https://flow.byu.edu/CCBlade.jl/stable/reference/#Output-Struct)

#Plots and Data analysis

#Plot Forces
figure()
plot(r/Rtip, out.Np)    #plots location (as percent of radius) and the normal force per unit length 
plot(r/Rtip, out.Tp) #plots "" and the tangential force per unit length 
xlabel("r/Rtip")
ylabel("distributed loads (N/m)")
legend(["flapwise (normal force", "lead-lag (tangential force)"])

#Plot induced velocities, can be used for future interaction calculations
figure()
plot(r/Rtip, out.u/Vinf) #axial induced velocity/Inflow velocity, so it is a percent, ie normalized
plot(r/Rtip, out.v/Vinf) #tangential induced velocity (swirl)/""...
xlabel("r/Rtip")
ylabel("Induced velocity at rotor disk (normalized)")
legend(["axial velocity", "swirl velocity"])


#Prep for Coefficient analysis (thrust, power efficiency etc) at different advance ratios
#Advance ratio = ratio of forward speed/ (rotational speed times prop diameter)

#Inputs
nJ = 20     #number of advance ratios to evaluate
J = range(.1, .6, length = nJ) #creates a range of advance ratios from .1 to .6, what are normal vals?
Omega = 5400.0*pi/30  #rpms to rad/sec 
n = Omega/(2*pi)    # converts radians persecond to just rotations per second, same as rpm/60
D = 2 * Rtip    #Diameter of the prop is 2* the radius, do I ignore the hub?

eff = zeros(nJ)     #creates arrays for efficiency, coef of thrust, and coef of torque 
CT = zeros(nJ)
Cq = zeros(nJ)      #coef of torque is requried torque over theoretical required torque

#Coef Solver

for i = 1:nJ
    local Vinf = J[i] * D * n   #makes a local inflow veloc var at each advance ratio 

    local op = simple_op.(Vinf, Omega, r, rho)  #creates op pts at each blade section/location
    outputs = solve.(Ref(rotor), sections, op) #uses all data from above plus local op conditions
    T, Q = thrusttorque(rotor, sections, outputs)   #calcs T & Q at each sec w/given conditions
    eff[i], CT[i], CQ[i] = nondim(T, Q, Vinf, Omega, rho, rotor, "propeller")
    # calcs the coef of each panel under the given conditions
end

#Prep to compare to experimental data 
exp = [
0.113   0.0912   0.0381   0.271
0.145   0.0890   0.0386   0.335
0.174   0.0864   0.0389   0.387
0.200   0.0834   0.0389   0.429
0.233   0.0786   0.0387   0.474
0.260   0.0734   0.0378   0.505
0.291   0.0662   0.0360   0.536
0.316   0.0612   0.0347   0.557
0.346   0.0543   0.0323   0.580
0.375   0.0489   0.0305   0.603
0.401   0.0451   0.0291   0.620
0.432   0.0401   0.0272   0.635
0.466   0.0345   0.0250   0.644
0.493   0.0297   0.0229   0.640
0.519   0.0254   0.0210   0.630
0.548   0.0204   0.0188   0.595
0.581   0.0145   0.0162   0.520
]