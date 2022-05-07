# AeroDynamic Term Dictionary

## Airfoil/Airplane Terms
**Angle of Attack-** Angle between the freestream direction and the chord line, usually relative to horizontal.
**Airfoil Polar-** Polars are graphs comparing the coefficients of lift and drag (among other things) with the angle of attack
**Stall-** When the angle of attack is too large the the flow separates from the airfoil and causes an increase in drag and a decrease in lift. Side note- downward twist towards the wing tips allows control surface near the tip of a wing to still have control even if the main portion of the wing has stalled.
**Pitch-** The angle with respect to horizontal, angle of attack is with reference to freestream direction, so in horizontal flight, these can be the same.
**Twist-** The angle between the root twist and tip twist, ie the tip can be twisted more up than the root, they can be reported separately if there is a linear twist profile as then we can know the twist anywhere along the airfoil.(For airplanes, typically angle from horizontal)
**Chord-** The length of the "chord line" (line connecting the leading and trailing edge of an airfoil).
**Camber-** Also called the maximum camber, it is the max distance between the Chord line and camber line (a line in the middle of the upper and lower, or Pressure and suction, surfaces of the airfoil). Higher camber ie more asymmetry between upper and lower surfaces generally increases lift.
**Airfoil Thickness-** The maximum thickness of an airfoil, or the distance between the upper and lower (Pressure and Suction) surfaces.
**Thickness-to-Chord ratio (*t/c*)-** The ratio of thickness of an airfoil divided by the chord (chord length) of the airfoil. Typically between 8-14% for typical aircraft airfoils
**NACA 4-series-** Airfoils defined by equations for thickness and camber and parameterized by four digits. The first is max camber as a % of the chord, the second is the location of the max camber in tenths of the chord, and the last two (really one number) are the thickness-to-chord ratio, ie the thickness of the airfoil in percent chord. Sometimes still used in simpler geometries, no longer the main airfoil parameterization method.
**Propeller Efficiency-** ratio of useful power-out from thrust, divided by the power in from shaft assuming no losses from motor power out. Basically exactly what you would expect for efficiency.
$$ \eta = \frac{P_{out}}{P_{in}} = \frac{TV}{Q\Omega} $$
*Explanation-* eta = Power output over power input, or Thrust times Velocity over Torque times rotational speed

**Advance Ratio (*J*)-** Ratio of forward speed to rotation speed times a constant
$$ J = \frac{V}{nD}  $$
*Explanation-* J = Forward velocity divided by n = some constant, D = rotation speed

**Tip Speed Ratio-** Ratio between the tip speed of the blade and the wind speed 
$$ \lambda = \frac{\omega R}{v}  $$
*Explanation-* Lamda = omega (rotational velocity) times by Propeller Radius divided by v = wind velocity (freestream)

**Freestream Velocity-** The velocity of the fluid around an airfoil before it interacts with the boundary layer of the airfoil, ie before it has been slowed down at all.
**Drag-** The force parallel to the freestream direction, resisting forward motion on the airfoil. Made up of Pressure drag and shear drag
**Lift-** The force perpendicular to the freestream direction that "lifts" the airfoil up.
**Flow Angle-** For rotating airfoils, it is the angle between the resultant wind speed (from wind and rotation) to the plan of rotation, also the sum of the angle of attack and the pitch angle


## Useful Coefficients
**Coefficient of Drag (*Cd*)-** It represents the drag effects of a certain fluid on an object. A lower coefficient means lower drag effects.
$$ C_d = \frac{F_d}{\frac{1}{2} \rho v^2A}  $$
*Explanation-* Cd = Coefficient of Drag, Fd is the force of drag, rho is the density of the fluid, v is the velocity of the fluid, A is the relevant surface area. 

**Coefficient of Lift (*Cl*)-** Theoretically 2pi, typically less in actuality. 
$$ C_l= m (\alpha - \alpha_0) $$
*Explanation-* m = lift curve slope from the linear portion of cl vs angle of attack graph.  alpha- angle of attack, alphanot- zero-lift angle of attack 

**Coefficient of Power (*Cp*)-**  It represents efficiency and is Power output of the turbine divided by the  total power from the fluid flowing into the blades. 
$$ C_p = \frac{P}{\rho n^3 D^5} $$
*Explanation-* P = Power, rho = density of the fluid, n = propeller rotation rate in revolutions per second, D = propeller diameter

**Coefficient of Torque (*CQ*)-**  Ratio of required torque to theoretical propeller required torque.
$$ C_Q = \frac {Q}{\rho n^2 D^5} $$
*Explanation-* Q = Torque, rho = density of the fluid, n = propeller rotation rate in revolutions per second, D = propeller diameter

**Coefficient of Thrust (*CT*)-**  The ratio of output thrust to propeller theoretical potential thrust.
$$ C_T = \frac {T}{\rho n^2 D^4} $$
*Explanation-* T = Thrust, rho = density of the fluid, n = propeller rotation rate in revolutions per second, D = propeller diameter

## Useful Dimensionless Numbers 
**Reynolds Number (*Re*)-** The Ratio of inertial force of the fluid over the viscous force of the fluid. For commercial airplanes the numbers will be in the millions to tens of millions (ie inertial forces play a larger part than viscous forces), for a UAV it will be in the 10,000 range.
$$ Re = \frac {\rho V c}{\mu} $$
*Explanation-* Rho = density of the fluid, V = velocity of the fluid (free stream), c = relevant length scale (the chord for an airfoil or aerodynamic chord for a  wing), mu = dynamic viscosity of the fluid

**Mach Number (*M*)-** Ratio between the freestream velocity over the speed of sound, ie when M = 1 = speed of sound, above = supersonic, above 5 = hyper sonic, below 1 = subsonic, below .3 = most likely incompressible flow and the Mach Number doesn't matter. 
$$ M = \frac{V}{a} $$
*Explanation-* V = freestream velocity, a = the speed of sound

## Github
**Repository-** Commonly called a "Repo", it is the location where all of the relevant files for a project are stored.
**Branch-** A branch is way to "branch off" of the main project without actually changing it, ie creating a new version that you can edit and then merge later
**Merge-** When you recombine a branch of a repo with the main repo
**Pull Request-** A pull request is a proposal to pull your branch or changes into another/the main branch. It allows you to view, comment, and discuss differences in the code of your branch vs the main branch.

## Other Terms
**Validation-** "Are we building the right product?" - Barry Boehm, comes after verification, it is looking to see if the outcomes meet the customer requirements and expected outcomes.
**Verification-** "Are we building the product right?" - Barry Boehm comes before validation and is to verify that the process in the code is correct and that there are no bugs, it is not actually running/executing the code, but instead is like "walking through" it and and checking the process
**Error-** Also called absolute error, this is the measure of how far off from the true/correct value you are. If the actual value is not known, it can be the uncertainty ie plus or minus a certain amount in the same units as the measurement
**Relative Error-** Shows what percent of the real value the absolute error is, ie how far off are we? So if the real value is 1 and we are off by .1 (absolute error) then relative error is .1/1 = 10%

## Sources
Dr. Ning's 415 (Flight Vehicle Design) Textbook
Several Google searches and random websites