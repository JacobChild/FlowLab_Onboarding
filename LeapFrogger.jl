#=
LeapFrogger.jl
Jacob Child
Flow Lab Onboarding Project
April 23rd, 2022

Equations
V(vector) = (Gamma(vector) cross r(vector))/(2*pi*r^2)

Pseudocode- Using the reference image and the right hand rule, compute the velocities for each 
vortex "point" individually. The bottom left point is 1, then go clockwise. 
For Pt1- compute the effects of pts 2,3,4 on it, add them together for Pt1TotVeloc, then using a 
small time step, compute its new location and repeat the process :||
All points have the same process
Things to watch- don't update the position of any points before all the calculations for
that time step are done.

Coordinate system:
(0; 0; 0) = [x; y; z] = inbetween and in line with pts 1 and 2
Pt 1 = btm left, the rest follow clockwise. 
Righthand rule = out of page = positive, into page = negative.
=#

#needed packages
using LinearAlgebra
using Plots
#using DataFrame

#Inputs/Variables
d = 1 #distance between vortex ring top and bottom points, and between rings
G = 1   #strength of vorticity = Gamma, same for all rings/points
dt = .01 #time step size
x = collect(0:dt:480) #gives a time vector from 0-xseconds with .01 second step size
plot()  #to clear the plot

#p1 = zeros(lenght(x), 3)
#Point 1
Pt1Gamma = [0; 0; -G]
Pt1Location = [0; -d/2; 0]
x1 = []
y1 = []

#Point 2
Pt2Gamma = [0; 0; G]
Pt2Location = [0; d/2; 0]
x2 = []
y2 = []

#Point 3
Pt3Gamma = [0; 0; G]
Pt3Location = [d; d/2; 0]
x3 = []
y3 = []

#Point 4
Pt4Gamma = [0; 0; -G]
Pt4Location = [d; -d/2; 0]
x4 = []
y4 = []


#before putting everything in functions, practice for pt 1
#2 on 1
# R21 = Pt1Location - Pt2Location #make sure positive vs negative is correct!
# println(R21)
# MagnitudeR21 = sqrt(R21[1]^2 + R21[2]^2 + R21[3]^2) #the magnitude of the distance btwn the 2 pts
# V21Inf = cross(Pt2Gamma, R21) / (2*pi*MagnitudeR21) #The influence pt2 has on 1 for velocity
# println(V21Inf)
# dt = .01    #time step
# Pt1dp = V21Inf * dt 
# Pt1Location = Pt1Location + Pt1dp
# println(Pt1Location)

function v_influence(Location1, Location2, Gamma)
    R = Location1 - Location2 #how to make sure this is right positive vs negative?
    MagnitudeR = sqrt(R[1]^2 + R[2]^2 + R[3]^2) #the magnitude of the distance btwn the 2 pts
    VInf = cross(Gamma, R) / (2*pi*MagnitudeR^2) #The influence pt2 has on 1 for velocity
    return VInf
end

function total_influence(VAInf,VBInf, VCInf)
    VTotInf = VAInf + VBInf +VCInf
    return VTotInf

end

function update_location(Location, VInfTot, dt)
    dp = VInfTot*dt
    NewLocation = Location + dp
    return NewLocation
end

###Looping Function

for i in x #i=1:length(x)

    #Point 1
    V21Inf = v_influence(Pt1Location, Pt2Location, Pt2Gamma)
    V31Inf = v_influence(Pt1Location, Pt3Location, Pt3Gamma)
    V41Inf = v_influence(Pt1Location, Pt4Location, Pt4Gamma)
    V1InfTot = total_influence(V21Inf, V31Inf, V41Inf)
    #println(V21Inf)

    #Point 2
    V12Inf = v_influence(Pt2Location, Pt1Location, Pt1Gamma)
    V32Inf = v_influence(Pt2Location, Pt3Location, Pt3Gamma)
    V42Inf = v_influence(Pt2Location, Pt4Location, Pt4Gamma)
    V2InfTot = total_influence(V12Inf, V32Inf, V42Inf)

    #Point 3
    V13Inf = v_influence(Pt3Location, Pt1Location, Pt1Gamma)
    V23Inf = v_influence(Pt3Location, Pt2Location, Pt2Gamma)
    V43Inf = v_influence(Pt3Location, Pt4Location, Pt4Gamma)
    V3InfTot = total_influence(V13Inf, V23Inf, V43Inf)

    #Point 4
    V14Inf = v_influence(Pt4Location, Pt1Location, Pt1Gamma)
    V24Inf = v_influence(Pt4Location, Pt2Location, Pt2Gamma)
    V34Inf = v_influence(Pt4Location, Pt3Location, Pt3Gamma)
    V4InfTot = total_influence(V14Inf, V24Inf, V34Inf)

    #Update Location of All points
    global Pt1Location = update_location(Pt1Location, V1InfTot, dt) #Make sure to update the positions last!
    # p1[i,:] = update_location(Pt1Location, V1InfTot, dt)
    global Pt2Location = update_location(Pt2Location, V2InfTot, dt)
    global Pt3Location = update_location(Pt3Location, V3InfTot, dt)
    global Pt4Location = update_location(Pt4Location, V4InfTot, dt)
    global x1 = push!(x1,Pt1Location[1]) #creates a series of x and y values for each point 
    global y1 = push!(y1,Pt1Location[2]) 
    global x2 = push!(x2,Pt2Location[1]) 
    global y2 = push!(y2,Pt2Location[2])
    global x3 = push!(x3,Pt3Location[1]) 
    global y3 = push!(y3,Pt3Location[2]) 
    global x4 = push!(x4,Pt4Location[1]) 
    global y4 = push!(y4,Pt4Location[2])
    #x1,y1,_ = update_location(Pt1Location, V1InfTot, dt)
    # #debugging
    # A = V12Inf[1]   #x direction velocity influence of pt 3 on pt 2
    # #B = V12Inf[2]   #y direction velocity influence of pt 3 on pt2
    # display(plot!([A], legend = false, m = 3))


end
#plots
#Vortex Ring 1 (back)
display(plot!(x1,y1,  label = false, color = :blue, xlabel="Position (m)", ylabel = "Position (m)")) #plot only does vectors, not floats etc, number values are plot commands
display(plot!(x2,y2,  label = false, color = :blue, title = "Leap Frogging Vortex Ring Positions Over Time"))

#Vortex Ring 2 (front)
display(plot!(x3,y3,  label = false, color = :green))
display(plot!(x4,y4,  label = false, color = :green))

# plt = plot(xaxis="X distance (m)", yaxis="Y distance (m)")
# plot!(x1, y1)
# display(plt)


#println(Pt1Location, Pt2Location, Pt3Location, Pt4Location)