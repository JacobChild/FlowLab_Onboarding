#=
LeapFrogger.jl
Jacob Child
Flow Lab Onboarding Project
April 23rd, 2022

Equations
V(vector) = (Gamma(vector) cross r(vector))/(2*pi*r^2)

=#

#test stuffs
Gamma = [0 0 1]
d = 1
t = [0:.01:.05;]
transpose(t)
r =  1   #not sure what the r vector is?

V = ((Gamma.^t)*r)/(2*pi*d^2/4)
#println(t)
println(V)

