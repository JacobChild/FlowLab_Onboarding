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
    println("\n t:", t)

    return c, p, t
end

"""
thickness_calculations(x,p)
Takes in the x coordinates, and max thickness (p), calculates thickness, and ouputs the array

"""
function thickness_calculations(x,p)
    m = p
    TArray = 10*m.*(.2969*sqrt.(x) - .1260*x - .3537*x.^2 + .2843*x.^3 - .1015*x.^4)
    return TArray
end

"""
zbar_calculations(c,p,x) 
Takes in c, p, and x values and outputs zbar 
"""
function zbar_calculations(c,p,x)
    for i in x
        if i <= p
            global zbar = c.*(2*p*x .- x.^2)/(p^2)
        elseif i > p
            global zbar = c.*(1 - 2*p .+ 2*p*x .- x.^2)/((1-p)^2)
        else 
            global zbar = "you wenti"
        end
    end

    return zbar 

end

#quick test function to understand iterating
function test_iterator(x)
    for i in x
        if i <= .3
            println(i+7)
        elseif i > .3
            println(i)
        end

    end
    
end

#quick test for zbar function
x = range(0,1, step=.1)
c =.2
p = .04
zbar = zbar_calculations(c, p, x)
println(zbar)
