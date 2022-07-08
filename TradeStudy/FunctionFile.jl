#FunctionFile.jl

#Remember for surfaces to work, it must return a matrix with rows and cols the size of the inputs

function RVarSurface(SRPM,NewRTip;spitout)

    for i = 1:NumofRadiis, j = 1:NumofVars
        if rvar[1].*NewRTip[i] < Rhub
            NewRHub = rvar[1].*NewRTip[i] #This is the new hub radius
            #println("Rhub was updated to $NewRHub for RVar")
        else
            NewRHub = Rhub
        end
        SOmega = (pi/30)*SRPM[j]  #converts rpm to rad/s, derivation: rpm*pi*360deg/(60sec*180deg)-> rpm*pi/30
        RVarRotor = Rotor(NewRHub, NewRTip[i], NumofBlades) #This is the new rotor object
        sections = Section.(rvar.*NewRTip[i], chordrvar.*NewRTip[i], twist, airfoils)
        op = simple_op.(Vinf, SOmega, rvar.*NewRTip[i], rho) 
        out = solve.(Ref(RVarRotor), sections, op)
        T, Q = thrusttorque(RVarRotor, sections, out)   #calcs T & Q at each sec w/given conds, sums them for the whole rotor
        EffRVar[i,j], CTRVar[i,j], CQRVar[i,j] = nondim(T, Q, Vinf, Omega, rho, RVarRotor, "helicopter")
    end
    if spitout == "Eff"
        println(EffRVar)
        return EffRVar
    elseif spitout == "CT"
        return CTRVar
    elseif spitout == "CQ"
        return CQRVar
    else
        println("Keyword argument error")
    end
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