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

    Goal = clean code, lots of well commented functions
=#


"""
user_input Function to take in the 4 NACA numbers
"""
function user_input()
    println("Input the 4 NACA numbers ie 2412 \nNACA ")
    input1 = readline()
    c = parse(Int64, input1[1])
    c = c//100  #converts the chord to the proper units
    p = parse(Int64, input1[2])
    p = p/10    #converts the chord position to the proper units
    t = parse(Int64, input1[3:4])
    t = t/100
    println("\n t:", t)

    return c, p, t, input1
end

"""
thickness_calculations(x,p)
Takes in the x coordinates, and max thickness (p), calculates thickness, and ouputs the array

"""
function thickness_calculations(x,t)
    m = t
    TArray = 10*m*(.2969*sqrt(x) - .1260*x - .3537*x^2 + .2843*x^3 - .1015*x^4)
    return TArray
end

"""
zbar_calculations(c,p,x) 
Takes in c, p, and x values and outputs zbar 
"""
function zbar_calculations(c,p,x)
    zbar = Float64[] #makes zbar an array of zeros the length of x, hopefully?
    for i in x
        if i <= p
             zbarInput = c*(2*p*i - i^2)/(p^2)
             push!(zbar,zbarInput)
            
        elseif i > p
             zbarInput = c*(1 - 2*p + 2*p*i - i^2)/((1-p)^2)
             push!(zbar,zbarInput)
             
        else 
             zbar = "you wenti"
             break
            
        end
        
    end
    return zbar
end

"""
zupper_lower(zbar, TArray)
calculates zu and zl and returns both
"""
function zupper_lower(zbar, TArray)

    zu = zbar .+ TArray/2
    zl = zbar .- TArray/2
    
    return zu, zl
end


#Actual code- See how it goes 

c, p, t, fromuser = user_input()  #function to get user input
#x = range(0,1,step=.1)  #creates an x range
x = collect(0:0.010:1)
TArray = thickness_calculations.(x, t)   #calculates the thicknesses
zbarn = zbar_calculations(c, p, x) #calculate and returns zbar
zu, zl = zupper_lower(zbarn, TArray) #calculates and returns zupper and z lower

#outputs

#println("\n \n", zu, "\n \n", zl) #for vs code

using DataFrames #only seems to work in the Julia REPL

DataFrame()

Data = DataFrame("X" => x, "Zupper"=> zu, "Zlower" => zl, "Thickness" => TArray, "Zbar" => zbarn)

display(Data)

#look into aspect ratio for plots 
 using Plots
 plt = plot(x,zu)   #saving the plot to an element
 plot!(x,zl)    #adding/updating to that plot
 #Plot settings
 plot!(xlims = (-.1, 1.1), title = ("NACA", fromuser), xlabel = "Position", ylabel = "Position")
 plot!(aspectratio = :equal)

 display(plt)   #displaying the plot by calling the saved element