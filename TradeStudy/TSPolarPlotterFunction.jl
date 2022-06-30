#=
TSPolarPlotterFunction.jl
Jacob Child
May 20th, 2022

High Level Project Overview-
Files in Project: ExtenderNSmoother.jl, PolarPlotter.jl, AirFoilCoordinatesFunction.jl 
PolarPlotter will be created and made to run as a stand alone file ("default" mode discussed below).
It will then be changed to call AirFoilCoordinatesFunction if need be, and AirFoilCoordinatesFunction
will be turned into a function rather than a stand alone file. ExtenderNSmoother will then be made to
run as it calls PolarPlotter, now turned into a function, for data.

Pseudocode-
Purpose- This code will generate the polar plots of the given Airfoil. It will ultimatley
probably export those results to the ExtenderNSmoother

Layout-
The user will choose whether to run from a loaded coordinates file, plot a new one, or use the default.
This can be done through state code with 3 states. Once the coordinates have been imported, XFoil will
be called and used  to generate the polar plots from -10deg to 15deg.

=#
function TSPolarPlotterFunction(x, y, Re, mach, ncrit, alpha) 

    #Solve for the Polars
    Xfoil.set_coordinates(x,y)
    Xfoil.pane(;npan = 160)    #this "panes" the coordinates, not sure what it does, smooths I think?
    
    #OUTPUTS
    n = length(alpha)
    Cl = zeros(n)
    Cd = zeros(n)
    Cdp = zeros(n) 
    Cm = zeros(n)
    Converged = zeros(Bool, n)
    
    #Calculations 
    #=
    for i = 1:n
        Cl[i], Cd[i], Cdp[i], Cm[i], Converged[i] = Xfoil.solve_alpha(alpha[i], Re; mach, iter = 100, ncrit, percussive_maintenance=true)
        #The above calculates the coefficients at each angle (i of alpha)
        #for some reason it isn't converging on the last one?
    end
    =#
    Cl, Cd, Cdp, Cm, Converged = Xfoil.alpha_sweep(x, y, alpha, Re; mach, iter = 100, ncrit, percussive_maintenance=false)
    @show Converged
    #Present Data TURNED OFF FOR FUNCTION 
     DataFrame()
    
    Data = DataFrame("Angle" => alpha, "Cl"=> Cl, "Cd" => Cd, "Cdp" => Cdp, "Cm" => Cm, "Convergence" => Converged)
    
    display(Data)
    
    return alpha, Cl, Cd, Cdp, Cm
    end