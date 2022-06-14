# CCBlade and BEM
**The Theory**
CCBlade runs off of the Blade Element Momentum theory. This theory basically says that if you draw a control volume around the rotating blade disk and do a momentum balance of the flow going in and out, that that momentum will be conserved. While that is true, there are several helpful assumptions made along the way that make the theoretical model tend to overpredict the capabilities of the actual propeller (see my report for more discussion)
**The How**
CCBlade is able to do this computation by discretizing along a propeller blade, defining the flow and blade characteristics at that section and computing its effects. It then repeats the process all along the propeller, sums the results and multiplies that by the number of blades. This makes CCBlade a relatively fast approximator of blade performance, and while not entirely accurate it provides excellent estimations that help in propeller analysis.

**The Capabilities**
CCBlade is incredibly powerful and can operate on a relatively simple level, or can become more complicated if wanted. At its simple level it is provided with the specifications of the rotor (diameter, twist, and airfoil cross section), the conditions of the incoming flow, and then computes at the specified number of discretized sections. At a higher level CCBlade is able to be told that a specific blade is made up of several different airfoil cross sections, where those cross sections occur, and what the inflow conditions are, and it is still able to compute the requested data. It is very quick and easy to analyze at a number of different RPMs and advance ratios and it overall helps to quickly analyze propellers.

**The Appendix and definitions:**
Of most use to me was understanding what was actually going in to and out of the CCBlade package function calls, so they will be explained below.

- Rotor: This is the first function call that comes up, it is "fed" The rotor info, so Rtip (radius of the blade to the tip, just the diameter of the blade/2, ie defined from the center of the hub/rotation), Rhub(radius of the blade hub), and B (number of blades). It returns the following as a stored object - (:Rhub, :Rtip, :B, :precone, :turbine, :mach, :re, :rotation, :tip), with only the 3 inputs mentioned above, precone = 0, turbine = false, mach, re, and rotation have nothing in them, and tip is PrandtlTipHub (*I am unsure of what this means, but maybe it is the type of tip or configuration?*). It will keep the rotor object to be called on for later use
- propgeom, while this is not a function, prop geometry is needed in the following formats: the various radius locations of the sections (ie take the column from 0 to 1 of the original prop file and multiply by Rtip), more = higher accuracy. It then needs the chord length at each radial location. It then needs the twist angle of the blade at each radial location (note, this must be in radians and UIUC is in degrees, so multiply by pi/180). *Note about the radial locations and chord:* If this data is downloaded off of the UIUC database, these numbers are all given in percentage of Radius, ie radial location / The radius and Chord / the radius, so to get it in the numbers CCBlade wants, the UIUC values will need to be multiplied by the radius.
- The next function call is AlphaAF. It is fed a file of airfoil data (made in xfoil) of the angles of attack and lift coefficients. The file format Is a title/header, the Reynolds number, and the mach number, each on their own row. That is followed by 3 columns of data , angle of attack, coef of lift, and coef of drag. Note, the angles of attack should be in radians, if not their is a keyword argument that can be changed. The data is then fit with an akima spline (to make it calculatable at any point) and is than stored all together as an object. AlphaAF does not have to take in a file, but can also take in each input by itself, ie AlphaAF(alpha, cl, cd, info, Re, Mach)
- Section this function defines each section at a radial location. It takes in the radial location, the chord length, the twist, and the aifoil data (from AlphaAF) for that section. If the blade has a constant airfoil crosssection you can define all sections at once by broadcasting with a "." and turning your airfoil data (for this example called af) into a refererence value, so 
```
sections = Sections.(r, chord, theta, Ref(af))
```

- Simple_op this defines the "operating points" for each radial location. An "operating point" is the conditions under which a certain radial section is at. This returns an object containing the operating point, and as above, if all can be set at once (ie all sections are under same conditions) the function can be broadcast with "." It takes in Simple_op(Vinf, Omega, r, rho;) then it has additional arguments for pitch, air viscosity (dynamic?), the speed of sound, and precone, but those all have default values unless changed. Vinf is the inflow velocity in the "x" direction. For a propeller The X direction is perpendicular to the rotating disk/rotation plane. Omega is the rpm in rad/s, r is the radial location, rho is the density of the air. *Note:* Simple_op operates under several assumptions; that Vx is just the freestream velocity/inflow velocity, that Vy is just \( \Omega * r \) 
(ie just the linear velocity of the rotor), this means there is no air inflow from the side, ie anything not perpendicular to the rotation disc. If these assumptions are not true, use OperatingPoint() this allows the user to input Vx and Vy etc...
- solve as the name of the function says, this is the solver, it takes in the rotor, section, and op structs made above and returns a struct of outputs, there is a large list, but the most notable are probably the normal force and tangential loading (VarName.Np and VarName.Tp respectively) as well as coefficients and several other things. it also can be "." broadcast to solve all the sections at once if it is the same, or several sections at once if there are multple airfoil section types and it needs to be run several times. note that each of the output structs is a vector of values, the length of the number of sections, for example out.Np[1] is the normal load for the first section of the blade etc...
- thrusttorque as the name implies, this calculates the Thrust and Torque, that is done by integrating the necessary values down the length of the blade. It takes in the rotor object, sections object, and output struct. It returns two float64s in the order of T,Q (Thrust, then Torque).
- nondim This nondimensionalizes the outputs (ie returns coefficients w/o units) the outputs depend on the input argument of the rotor type (ie windturbine vs propeller vs helicopter) for a propeller it returns efficiency, Coef of Thrust, and Coef of torque (CQ) in that order. It takes in the integrated thrust, torque, Vhub (a hub speed for turbine normalization in m/s, this is just Vinf), Omega (the rpms in rad/s), rho (density), the rotor object, and the string rotortype (ie "propeller"). Note: by using a for loop you can calculate these values depending on an array of Advance Ratios (J)

### Definitions
**Coefficient Definitions in my own words**
*Specifically for when plotted against the Advance Ratio*
- **Advance Ratio**
$$ J = \frac{V}{nD}  $$
*Explanation-* J = Forward velocity divided by n = revolutions per second, D = Diameter
*Definition:* The Advance Ratio is the ratio between the forward velocity and the velocity of the rotating blade tip. Basically how fast is it moving forward in relation to how fast the blade is spinning (ie an increasing advance ratio would be moving forward faster compared to how fast the rpms are, you can imagine an rc plane speeding up with a fixed rpm)
- **Coefficient of Thrust**
- $$ C_T = \frac {T}{\rho n^2 D^4} $$
*Explanation-* T = Thrust, rho = density of the fluid, n = propeller rotation rate in revolutions per second, D = propeller diameter
*Definition:* The coefficient of thrust is how much thrust is produced over the theoretical amount of thrust that could be produced. We know that the faster a propeller spins the more thrust it produces, however, when plotted against increasing advance ratios, the Coef of Thrust goes down, this is because *in relation to* the air flowing in, the propeller produces less thrust. Thrust comes from grabbing air and throwing it back faster, when the air comes in fast, and you only have a fixed speed you can throw it back at, the thrust that can be produced goes down. The thrust potential of the blade (the bottom of the equation) stays the same, but the Thrust produced decreases, so it goes down.
- **Coefficient of Torque**
$$ C_Q = \frac {Q}{\rho n^2 D^5} $$
*Explanation-* Q = Torque, rho = density of the fluid, n = propeller rotation rate in revolutions per second, D = propeller diameter
*Definition:* The amount of torque required to turn the blade, normalized in a way that makes it non dimensional, useful in that it can be converted to Coefficient of Power
- **Derivation of Coefficient of Torque to Coefficient of Power**

$$C_Q = \frac {Q}{\rho n^2 D^5} $$
$$ P = Q\omega $$
$$ C_p = \frac{P}{\rho n^3 D^5} \rArr \frac{Q\omega}{\rho n^3 D^5} $$
$$ n = \frac {\omega}{2\pi} $$
$$ C_p \rArr \frac {Q \omega 2\pi}{\rho \omega n^2 D^5} \rArr  \frac{2\pi Q}{\rho n^2 D^5} = 2\pi C_Q $$
$$  C_p = 2\pi C_Q $$


- **Coefficient of Power**
$$ C_p = \frac{P}{\rho n^3 D^5} $$
*Explanation-* P = Power, rho = density of the fluid, n = propeller rotation rate in revolutions per second, D = propeller diameter
*Definition:* The Coefficient of Power is the Power that is produced by the propeller divided by the theoretical amount of power that could be produced by that propeller at that rpm etc... When plotted against the advance ratio it peaks when the prop is both producing and drawing the most power, and goes down on both sides from there. 

- **Efficiency**
$$ \eta = \frac{P_{out}}{P_{in}} = \frac{TV}{Q\Omega} $$
*Explanation-* eta = Power output over power input, or Thrust times Velocity over Torque times rotational speed
*Definition:* Efficiency is simply how much power produced for the power coming in, ie 80% efficiency means 80% of the power coming in is getting converted to power out. As the advance ratio increases, Efficiency typically increases as the motor is working less hard to speed the air up. 
