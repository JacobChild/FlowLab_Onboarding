#OptimizationFunctionFile.jl

#Remember for surfaces to work, it must return a matrix with rows and cols the size of the inputs
#The output var names were changed to differ from the normal plots by putting an "S" at the end, 
#   ie FMRVarS etc, The "S" represents "Surface"

```
New Functions used for the optimization
```
function InitialPopulation(PopSize, NumofVars, lower, upper)
   #Make a vector of Individuals
   OutputPopInit = fill(zeros(NumofVars),PopSize) #makes a vector of vectors
	for ind in 1:PopSize
	Individual = zeros(NumofVars)
		for var in 1:NumofVars
		Individual[var] = rand(range(lower[var],upper[var],length=500))
		end
	OutputPopInit[ind] = Individual
	end
    return OutputPopInit

end

function ConstraintFunction(x)
    Q = zeros(1) #this makes it a 1 element vector
    T = zeros(1)
    #println("Constraint Input $x")
    CurrentChordDist = x[1:round(Int64,NumofVars/2)]
    CurrentTwistDist = x[round(Int64,NumofVars/2)+1:end]
    sections = Section.(r, CurrentChordDist, CurrentTwistDist, airfoils)
    op = simple_op.(Vinf, Omega, r, rho)
    out = solve.(Ref(MyRotor), sections, op)
    T[1], Q[1] = thrusttorque(MyRotor, sections, out)   #calcs T & Q at each sec w/given conds, sums them for the whole rotor

    return Q #Returns the Torque 
end

function ObjectiveFunction(x)
    
    CurrentChordDist = x[1:round(Int64,NumofVars/2)]
    CurrentTwistDist = x[round(Int64,NumofVars/2)+1:end]
    sections = Section.(r, CurrentChordDist, CurrentTwistDist, airfoils)
    op = simple_op.(Vinf, Omega, r, rho)
    out = solve.(Ref(MyRotor), sections, op)
    T, Q = thrusttorque(MyRotor, sections, out)   #calcs T & Q at each sec w/given conds, sums them for the whole rotor
        
    OutputT = -1.0*T
    return OutputT #Returns the Thrust, negative so it can "minimize" the function 
end

function Verification(x)
    CurrentChordDist = x[1:round(Int64,NumofVars/2)]
    CurrentTwistDist = x[round(Int64,NumofVars/2)+1:end]
    sections = Section.(r, CurrentChordDist, CurrentTwistDist, airfoils)
    op = simple_op.(Vinf, Omega, r, rho)
    out = solve.(Ref(MyRotor), sections, op)
    T, Q = thrusttorque(MyRotor, sections, out)   #calcs T & Q at each sec w/given conds, sums them for the whole rotor
    OutputT = -1.0*T
    #println("Thrust = $T")
    #println("Torque = $Q")
    return T, Q, OutputT
end

function WhatChanged(New,chord,twist)
    #Chord Changes
    FPercentChordChange = (New[1:27] .- chord) ./ chord .* 100
    ActualChordChangeMM = (New[1:27] .- chord) .* 1000
    
    #Twist Changes
    FPercentTwistChange = (New[28:54] .- twist) ./ twist .* 100
    NewDegs = New[28:54].* 180 ./ pi
    OldDegs = twist .* 180 ./ pi
    ActualTwistChangeDeg =NewDegs - OldDegs
    
    #Plot
    PercentChangesPlot = plot(FPercentChordChange, margin = 17mm, xlabel = "Blade section", ylabel = "Percent Change", label = "Chord Change", legend =:right, title = "Percent Blade Geometry Change", grid=true, minorgrid=true)
    plot!(twinx(), FPercentTwistChange, linecolor =:orange, ylabel = "Percent Change", label = "Twist change", legend =:topright)
    
    ActualChangesPlot = plot(ActualChordChangeMM, margin = 17mm, xlabel = "Blade section", ylabel = "Chord Change (MM)", label = "Chord Change", legend =:right, title = "Blade Geometry Change", grid=true, minorgrid=true)
    plot!(twinx(), ActualTwistChangeDeg, linecolor =:orange, ylabel = "Twist Change (Deg)", label = "Twist change", legend =:topright)
    
    ActualGeometryPlot = plot(New[1:27] .* 100, margin = 17mm, xlabel = "Blade section", ylabel = "Chord (CM)", label = "Chord", legend =:right, title = "Blade Geometry", grid=true, minorgrid=true)
    plot!(twinx(), New[28:54] .* 180 ./ pi, linecolor =:orange, ylabel = "Twist (Deg)", label = "Twist", legend =:topright)
    
    OriginalGeometryPlot = plot(chord .* 100, margin = 17mm, xlabel = "Blade section", ylabel = "Chord (CM)", label = "Chord", legend =:right, title = "Original Blade Geometry", grid=true, minorgrid=true)
    plot!(twinx(), twist .* 180 ./ pi, linecolor =:orange, ylabel = "Twist (Deg)", label = "Twist", legend =:topright)
    
    display(PercentChangesPlot)
    display(ActualChangesPlot)
    display(ActualGeometryPlot)
    display(OriginalGeometryPlot)
    
    end

```
Old Functions needed for trade study and older etc
```

function RVarSurface(SRPM,NewRTip;spitout)
    for k = 1:length(SRPM)

        for i = 1:NumofRadiis
            if rvar[1].*NewRTip[i] < Rhub
                NewRHub = rvar[1].*NewRTip[i] #This is the new hub radius
                #println("Rhub was updated to $NewRHub for RVar")
            else
                NewRHub = Rhub
                #println("Rhub stays the same")
            end
            SOmega = (pi/30)*SRPM[k]  #converts rpm to rad/s, derivation: rpm*pi*360deg/(60sec*180deg)-> rpm*pi/30
            RVarRotor = Rotor(NewRHub, NewRTip[i], NumofBlades) #This is the new rotor object
            sections = Section.(rvar.*NewRTip[i], chordrvar.*NewRTip[i], twist, airfoils)
            op = simple_op.(Vinf, SOmega, rvar.*NewRTip[i], rho) 
            out = solve.(Ref(RVarRotor), sections, op)
            T, Q = thrusttorque(RVarRotor, sections, out)   #calcs T & Q at each sec w/given conds, sums them for the whole rotor
            FMRVarS[i,k], CTRVarS[i,k], CQRVarS[i,k] = nondim(T, Q, Vinf, SOmega, rho, RVarRotor, "helicopter")
        end
    end
    if spitout == "FM"
        return FMRVarS
    elseif spitout == "CT"
        return CTRVarS
    elseif spitout == "CQ"
        return CQRVarS
    else
        println("Keyword argument error")
    end
end

function CVarSurface(SRPM, VarChordPercent; spitout)
    for k = 1:length(SRPM)
        for i = 1:NumofChords
            SOmega = (pi/30)*SRPM[k]  #converts rpm to rad/s, derivation: rpm*pi*360deg/(60sec*180deg)-> rpm*pi/30
            sections = Section.(r, ((CVarchord.+VarChordPercent[i]).*Rtip), twist, airfoils)
            op = simple_op.(Vinf, SOmega, r, rho)
            out = solve.(Ref(MyRotor), sections, op)
            T, Q = thrusttorque(MyRotor, sections, out)   #calcs T & Q at each sec w/given conds, sums them for the whole rotor
            FMCVarS[i,k], CTCVarS[i,k], CQCVarS[i,k] = nondim(T, Q, Vinf, SOmega, rho, MyRotor, "helicopter")
            # calcs the coef of the blade under the given conditions at each advance ratio
        end
    end
    if spitout == "FM"
        return FMCVarS
    elseif spitout == "CT"
        return CTCVarS
    elseif spitout == "CQ"
        return CQCVarS
    else
        println("Keyword argument error")
    end
end


function PVarSurface(SRPM, SVarPitch; spitout)
    for k = 1:NumofVars
        
        for i = 1:NumofPitches
            SOmega = (pi/30).*SRPM[k]
            pitch = SVarPitch[i]*pi/180 #This is the pitch angle to use
            op = simple_op.(Vinf, SOmega, r, rho; pitch)
            out = solve.(Ref(MyRotor), PVarSections, op)
            T, Q = thrusttorque(MyRotor, PVarSections, out)   #calcs T & Q at each sec w/given conds, sums them for the whole rotor
            FMPVarS[i,k], CTPVarS[i,k], CQPVarS[i,k] = nondim(T, Q, Vinf, SOmega, rho, MyRotor, "helicopter")
            # calcs the coef of the blade under the given conditions at each advance ratio
        end
    end
    if spitout == "FM"
        return FMPVarS
    elseif spitout == "CT"
        return CTPVarS
    elseif spitout == "CQ"
        return CQPVarS
    else
        println("Keyword argument error")
    end
end

function JVarSurface(SRPM, VarAdvanceRatio; spitout)
    for k = 1:NumofVars
        nS[k] = SRPM[k]/60 #ie n = rotations per second 
        for i = 1:NumofAdvanceRatios
            local Vinf = VarAdvanceRatio[i] * D * nS[k]   #makes a local inflow veloc var at each advance ratio 
            local op = simple_op.(Vinf, Omega, r, rho)  #creates op pts at each blade section/location
            outputs = solve.(Ref(MyRotor), Normalsections, op) #uses all data from above plus local op conditions
            T, Q = thrusttorque(MyRotor, Normalsections, outputs)   #calcs T & Q at each sec w/given conds, sums them for the whole rotor
            FMJVarS[i,k], CTJVarS[i,k], CQJVarS[i,k] = nondim(T, Q, Vinf, Omega, rho, MyRotor, "helicopter")
            # calcs the coef of the blade under the given conditions at each advance ratio
        end
    end

    if spitout == "FM"
        return FMJVarS
    elseif spitout == "CT"
        return CTJVarS
    elseif spitout == "CQ"
        return CQJVarS
    else
        println("Keyword argument error")
    end

end


function FindFunc(Search::StepRangeLen,Find)
	stepsize= parse(Float64, string(Search.step)[30:findall(k->k==',',string(Search.step))[1]-1])/2
	location = findall(k-> k-stepsize<Find<k+stepsize,Search)

	return location
end


function MyTestFunc2(a,b; spitout) #? For reference only, practicing nested for loops
	println(length(a))	
z = zeros(length(a),length(b))
	println(x)
	for i = 1:length(a),j = 1:length(y)
       	if spitout == "add"
           		z[i,j] = a[i]+b[j]
       	elseif spitout == "hyp"
           		z = sqrt.(a.^2+b.^2)
       	elseif spitout == "wiggles"
           		z = sin(a) + cos(b) + cos(a)*sin(b)
       	else
           		z = 0
       	end
    end
	println(z)
       return z

end

function MyTestFunc(a,b; spitout) #? For reference only, practicing nested for loops on one line
	println(length(a))	
z = zeros(length(a),length(b))
	println(x)
	
       	if spitout == "add"
           		z = [a+b for a in x, b in y]
       	elseif spitout == "hyp"
           		z = sqrt.(a.^2+b.^2)
       	elseif spitout == "wiggles"
           		z = sin(a) + cos(b) + cos(a)*sin(b)
       	else
           		z = 0
       	end

	println(z)
       return z

end