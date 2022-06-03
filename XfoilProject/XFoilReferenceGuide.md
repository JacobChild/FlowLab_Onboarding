# XFoil Reference Guide

XFoil is made for the analysis and design of airfoils and streamlined bodies.

**Navigation-**

a "." or typing a command/function takes you down a menu and pressing enter/return takes you up

Typing a "?" shows current menu options for the location you are at.

**General Airfoil Analysis Process-**
1. Load an airfoil, if a NACA airfoil, type NACA and the defining numbers `Naca 2412`
Returns- Info about the airfoil including thickness, number of points/panels, etc.
2. Define conditions/settings for airfoil analysis
	a. `oper` this takes you within the "operating points" menu to begin defining conditions. returns `.OPERi` the i means "inviscid flow"
	b. change to viscous flow by typing `v` and then enter a desired Reynolds number
	c. to change the Mach number type `m`. Other parameters such as Ncrit are under the `vpar` menu. type `vpar` and then `n` to change the Ncrit value
	d. To get out of the `vpar` menu and back up to the `.OPERv` menu, press enter/return
	e. Xfoil solves through iterations seeking to converge (see explanation in the appendix) its default iteration number is very low, type `iter` and change the value to above 100 (I typically do 500)
3. Calculate the Polars etc... type `alfa` followed by the angle of attack you want to analyze, or if you want to analyze a sequence type `aseq` followed by 3 numbers, starting angle, stopping angle, and step size, so `aseq 0 10 1` would be from 0 to 10 degrees analyzing at every 1 degree. This returns a plot of Cl and alpha, as well as a table of the various coefficient values with alpha
4. To see a plot of all the polars and compare different airfoils
	a. Turn on point accumulation (basically it will remember what it calculates) type `pacc` it will then prompt for a file to save to and a file to dump to, these are optional and can be skipped by pressing enter. 
	b. Rerun your analysis as in step 3.
	c. Type `pplo` to see the plotted polars. 
	d. To compare different airfoils, type `pacc` to stop accumulating points momentarily. and go up a menu by pressing enter. Type in your new airfoil `Naca 4412` and type `oper`
	e. turn accumulation back on `pacc` and reanalyze by doing step 3.
	f. type `pplo` to see the new polar plot, as well as the old saved one on the same plot.

## Appendix
The XFoil Documentation can be found [here](https://web.mit.edu/drela/Public/web/xfoil/xfoil_doc.txt) and the program itself can be downloaded etc from [here](https://web.mit.edu/drela/Public/web/xfoil/)

### Goals -
To both analyze and design airfoils. 
The analysis of the airfoils can be done under various user set conditions and can return information ranging from Polar Plots to pressure distribution graphs.
The design of the airfoils can be done through defining geometric parameters (ie designing the airfoil), or defining the flow and surface speed conditions and creating an airfoil that would generate the given flow characteristics (ie defining the flow).

### Capabilities -
Analysis Capabilities-
Viscous and Inviscid flow, "limited trailing edge separation, lift and drag predictions just beyond CLmax". Corrects for some compressibility effects. Can vary Re, M, Ncrit, etc. Can plot polars, coefficients, and pressure distributions as well as streamline approximations. Input airfoils can be custom and from a file.

Analysis Limitations-
 XFoil doesn't do well if there is flow seperation, ie non streamlined or high angle of attack. It kind of runs off of the vortex panel method, but has lots more math over my head.
![b54bfad530aff0376c9fd97fc9e787f7.png](:/ecccb18476754e9bacf5fbde319eb886)

It predicts well up to stall point, then in the non-linear regions it is not as good, as well not as good with high value reynolds numbers
It can only analyze simple lift devices, ie airfoils, airfoils with simple flaps, etc. Nothing super complicated

Design Capabilities-
Both direct geometry design and full and mixed-inverse flow design. Airfoil blending, Flap deflection force analysis (good for knowing needed servo strength) 

Design Limitations-
Not a CAD program etc...


### Definitions -
Panel Method- Discretizing an airfoil/wing by breaking it up into "panels". The boundary conditions are specified, like know flow through the surface, and freestream velocity at a certain distance etc, and then each panel is given a "fluid element" (typically a vortex) that will create those flow conditions. Essentially each panel will have a point at which there is a vortex creating the flow conditions desired, and a point at which calculations are done to verify the flow conditions (called the collocation point). The effects of all of the panels can then be summed to give the entire flow field around the surface with the desired conditions. It can be applied to both 2d and 3d objects and helps with complex calculations. It is important to remember that the panels must be sufficiently small to provide the level of accuracy you want. An awesome resource can be found [here](https://open.oregonstate.education/intermediate-fluid-mechanics/chapter/the-panel-method-an-introduction/)

Ncrit- An amplification factor for the smoothness of the airflow approaching the airfoil. 

Negative Stall- When analyzing an airfoil from -pi to pi, the airfoil does a full 360 with respect to the flow. A normal stall is when the airfoil stops generating lift at a high angle of attack. When it passes that point it begins to generate "negative lift" and is sucked down. As it continues to rotate it will experience negative stall where it goes from generating "negative lift" to positive lift again.

Viterna Method- An extrapolation method built into the CCBlade.jl wrapper. It is made for airfoils and extrapolates the Cl and Cd data from -pi to pi [documentation](https://flow.byu.edu/CCBlade.jl/stable/reference/)

1D Spline- This is what is used to smooth the given input data. It is found in the Dierckx julia package and creates a 1D spine to be overlayed the input data. It is made with user given order and smooth factor inputs.

Linear Region- The region of angles of attack in which the Coefficients of Lift and drag increase in a linear manner. This is pre stall.

Angle of Attack Sweep- A function `aseq` that allows you to "sweep" through several angles of attack and perform the analysis calculations for an airfoil.

Cdp- Coefficient of Pressure Drag = Coef of Drag - Coef of Friction Drag

Cm- Coefficient of Moment, The ratio of the moment to the Lift. The moment is caused by the center of pressure not being at the center of the airfoil. In Xfoil or M/q, conventional = moment/(q*c^2)
$$ Cm = 1/2 * \rho xv^2 $$


## Resources
[Playlist for an aerodynamics class, includes xfoil and xflyer tutorials](https://www.youtube.com/playlist?list=PLBcnfVKyTTtLGW7H1ofYaRX_fz9Gr6V9g)

