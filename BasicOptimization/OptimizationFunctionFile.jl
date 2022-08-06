#OptimizationFunctionFile.jl

#Remember for surfaces to work, it must return a matrix with rows and cols the size of the inputs
#The output var names were changed to differ from the normal plots by putting an "S" at the end, 
#   ie FMRVarS etc, The "S" represents "Surface"

```
New Functions used for the optimization
```
function ConstraintFunction(x)
    CurrentChord = x[1]
    pitch = x[2]*pi/180
    sections = Section.(r, ((CVarchord.+CurrentChord).*Rtip), twist, airfoils)
    op = simple_op.(Vinf, Omega, r, rho; pitch)
    out = solve.(Ref(MyRotor), sections, op)
    T, Q = thrusttorque(MyRotor, sections, out)   #calcs T & Q at each sec w/given conds, sums them for the whole rotor
    #FMCVar[i], CTCVar[i], CQCVar[i] = nondim(T, Q, Vinf, Omega, rho, MyRotor, "helicopter")
    # calcs the coef of the blade under the given conditions at each advance ratio
    return Q #Returns the Torque 
end

function ObjectiveFunction(x)
    CurrentChord = x[1]
    pitch = x[2]*pi/180
    sections = Section.(r, ((CVarchord.+CurrentChord).*Rtip), twist, airfoils)
    op = simple_op.(Vinf, Omega, r, rho; pitch)
    out = solve.(Ref(MyRotor), sections, op)
    T, Q = thrusttorque(MyRotor, sections, out)   #calcs T & Q at each sec w/given conds, sums them for the whole rotor
    #FMCVar[i], CTCVar[i], CQCVar[i] = nondim(T, Q, Vinf, Omega, rho, MyRotor, "helicopter")
    # calcs the coef of the blade under the given conditions at each advance ratio
    return -T #Returns the Thrust, negative so it can "minimize" the function 
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