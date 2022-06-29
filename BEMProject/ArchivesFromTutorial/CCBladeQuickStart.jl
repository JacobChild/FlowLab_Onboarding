#=
CCBladeQuickStart.jl
Jacob Child 
May 30, 2022

Following the tutorial - https://flow.byu.edu/CCBlade.jl/stable/tutorial/

PseudoCode- from UIUC do an analysis of a thin 10 x 5 propeller, 10 means diameter of 10inches

cd("C:/Users/child/Documents/Flow_Lab/Onboarding/BEMProject")
=#

#Packages
using CCBlade, PyPlot

#Defining input values
Rtip = 10/2.0 * 0.0254  # inches to meters
Rhub = 0.10*Rtip
B = 2  # number of blades

rotor = Rotor(Rhub, Rtip, B)

propgeom = [0.15   0.130   32.76
0.20   0.149   37.19
0.25   0.173   33.54
0.30   0.189   29.25
0.35   0.197   25.64
0.40   0.201   22.54
0.45   0.200   20.27
0.50   0.194   18.46
0.55   0.186   17.05
0.60   0.174   15.97
0.65   0.160   14.87
0.70   0.145   14.09
0.75   0.128   13.39
0.80   0.112   12.84
0.85   0.096   12.25
0.90   0.081   11.37
0.95   0.061   10.19
1.00   0.041   8.99
]

r = propgeom[:,1] * Rtip    #sets up the varius radius sections
chord = propgeom[:,2] * Rtip    #normalizes the chord, what does that mean
theta = propgeom[:,3] * pi/180  #converts to radians

af = AlphaAF("naca4412.dat")   #are we in the right directory?

sections = Section.(r, chord, theta, Ref(af))   #this populates all the sections at once due to the "." broadcasting
    #I think Ref(af) makes af a "reference value" or something so it doesn't change, ie it isn't an array and everything else is

#Define the operating points, this means setting up the inflow and conditions
Vinf = 5.0
Omega = 5400 * pi/30    #rpms into rads/sections
rho = 1.225 #density of air

op = simple_op.(Vinf, Omega, r, rho)

out = solve.(Ref(rotor), sections, op)

out.alpha*180/pi    #out is a structure, thus this is converting alpha to rads 

#Plot outputs

figure(1)
plot(r/Rtip, out.Np)
plot(r/Rtip, out.Tp)
xlabel("r/Rtip")
ylabel("distributed loads (N/m)")
legend(["flapwise", "lead-lag"])

figure(2)
plot(r/Rtip, out.u/Vinf)
plot(r/Rtip, out.v/Vinf)
xlabel("r/Rtip")
ylabel("(normalized) induced velocity at rotor disk")
legend(["axial velocity", "swirl velocity"])


nJ = 20  # number of advance ratios

J = range(0.1, 0.6, length=nJ)  # advance ratio

Omega = 5400.0*pi/30
n = Omega/(2*pi)
D = 2*Rtip

eff = zeros(nJ)
CT = zeros(nJ)
CQ = zeros(nJ)

for i = 1:nJ
    local Vinf = J[i] * D * n

    local op = simple_op.(Vinf, Omega, r, rho)
    outputs = solve.(Ref(rotor), sections, op)
    T, Q = thrusttorque(rotor, sections, outputs)
    eff[i], CT[i], CQ[i] = nondim(T, Q, Vinf, Omega, rho, rotor, "propeller")

end

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
Jexp = exp[:, 1]
CTexp = exp[:, 2]
CPexp = exp[:, 3]
etaexp = exp[:, 4]


figure()
plot(J, CT)
plot(J, CQ*2*pi)
plot(Jexp, CTexp, "ko")
plot(Jexp, CPexp, "ko")
xlabel(L"J")
legend([L"C_T", L"C_P", "experimental"])

figure()
plot(J, eff)
plot(Jexp, etaexp, "ko")
xlabel(L"J")
ylabel(L"\eta")
legend(["CCBlade", "experimental"])
