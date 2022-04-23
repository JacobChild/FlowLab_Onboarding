#=
AirFoilCoordinatesArraysActivity.jl
Jacob Child
April 23rd, 2022

Equations Used
Zu(Zupper) = zbar + t/2
Zl(Zlower) = zbar -t/2
t = 10m(.2969sqrt(x) - .1620x - .3537x^2 + .2843x^3 - .1015x^4)
if X<=p
        zbar = c*(2*p*x-x^2)/(p^2)
    elseif x > p 
        zbar = c*(1-2*p+2*p*x-x^2)/((1-p)^2)
    end

Pseudocode- Take in NACA 4- Series Airfoil numbers- ie max camber (c),
    max camber position (p), thickness (t), and an array (x) to plot at.
    Do that by asking the user for c,p,t, then creating x (hardcode it for simplicity)
    Use the functions previously written to calculate zu and zl.


=#