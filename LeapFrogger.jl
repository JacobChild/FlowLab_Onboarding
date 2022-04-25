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

#Inputs/Variables
d = 1 #distance between vortex ring top and bottom points, and between rings
G = 1   #strength of vorticity = Gamma, same for all rings/points

#Point 1
Pt1Gamma = [0; 0; -G]
Pt1Location = [0; -d/2; 0]
Pt1TotVeloc = [0; 0; 0]

#Point 2
Pt2Gamma = [0; 0; G]
Pt2Location = [0; d/2; 0]
Pt2TotVeloc = [0; 0; 0]

#Point 3
Pt3Gamma = [0; 0; G]
Pt3Location = [d; d/2; 0]
Pt3TotVeloc = [0; 0; 0]

#Point 4
Pt4Gamma = [0; 0; -G]
Pt4Location = [d; -d/2; 0]
Pt4TotVeloc = [0; 0; 0]
