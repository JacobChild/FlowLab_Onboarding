#FunctionFile.jl

#Remember for surfaces to work, it must return a matrix with rows and cols the size of the inputs

function RVarSurface(SRPM,NewRTip;spitout)
    for k = 1:length(SRPM)

        for i = 1:NumofRadiis
            if rvar[1].*NewRTip[i] < Rhub
                NewRHub = rvar[1].*NewRTip[i] #This is the new hub radius
                println("Rhub was updated to $NewRHub for RVar")
            else
                NewRHub = Rhub
            end
            SOmega = (pi/30)*SRPM[k]  #converts rpm to rad/s, derivation: rpm*pi*360deg/(60sec*180deg)-> rpm*pi/30
            RVarRotor = Rotor(NewRHub, NewRTip[i], NumofBlades) #This is the new rotor object
            sections = Section.(rvar.*NewRTip[i], chordrvar.*NewRTip[i], twist, airfoils)
            op = simple_op.(Vinf, SOmega, rvar.*NewRTip[i], rho) 
            out = solve.(Ref(RVarRotor), sections, op)
            T, Q = thrusttorque(RVarRotor, sections, out)   #calcs T & Q at each sec w/given conds, sums them for the whole rotor
            FMRVar[i,k], CTRVar[i,k], CQRVar[i,k] = nondim(T, Q, Vinf, SOmega, rho, RVarRotor, "helicopter")
        end
    end
    if spitout == "FM"
        return FMRVar
    elseif spitout == "CT"
        return CTRVar
    elseif spitout == "CQ"
        return CQRVar
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
            out = solve.(Ref(rotor), sections, op)
            T, Q = thrusttorque(rotor, sections, out)   #calcs T & Q at each sec w/given conds, sums them for the whole rotor
            FMCVar[i,k], CTCVar[i,k], CQCVar[i,k] = nondim(T, Q, Vinf, SOmega, rho, rotor, "helicopter")
            # calcs the coef of the blade under the given conditions at each advance ratio
        end
    end
    if spitout == "FM"
        return FMCVar
    elseif spitout == "CT"
        return CTCVar
    elseif spitout == "CQ"
        return CQCVar
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
            out = solve.(Ref(rotor), PVarSections, op)
            T, Q = thrusttorque(rotor, PVarSections, out)   #calcs T & Q at each sec w/given conds, sums them for the whole rotor
            FMPVar[i,k], CTPVar[i,k], CQPVar[i,k] = nondim(T, Q, Vinf, SOmega, rho, rotor, "helicopter")
            # calcs the coef of the blade under the given conditions at each advance ratio
        end
    end
    if spitout == "FM"
        return FMPVar
    elseif spitout == "CT"
        return CTPVar
    elseif spitout == "CQ"
        return CQPVar
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