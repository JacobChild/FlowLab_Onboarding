# NACA 2412 Data and Discussion

### Introduction

The NACA 2412 airfoil was created, plotted, and then analyzed under different conditions through various programs In order to better understand basic principles of airfoil analysis and gain a familiarity with XFoil and the XFoil Julia wrapper. The results from the analyses are given below followed by validation and verification findings and discussions and then discussions of the overall findings from the data.

### Table of Contents

- Findings and Plots
  
  - M = 0, Ncrit = 7, Varied Reynolds numbers
  
  - M = 0.2, Ncrit = 11, Varied Reynolds numbers

- Validation and Verification 
  
  - Validation methods and error discussion
  
  - Verification methods and error discussion

- Overall findings discussion
  
  - Reynolds Number trends
  
  - Angle of attack trends
  
  - Stall trends

- Conclusion

- Appendix
  
  - High level code explanation and overview
  
  - Source code

<img title="" src="file:///C:/Users/child/Documents/Flow_Lab/Onboarding/XfoilProject/OutputImages/NACA2412PlotFig.jpg" alt="NACA2412Plot" data-align="center" width="332">

*Figure 1: A plot of the NACA 2412 Airfoil. Generated using cosine spacing and 20 points on both the top and bottom (for the source code see [Here](https://github.com/JacobChild/FlowLab_Onboarding/blob/9f6389461f356a878bb87c64cd1cda8c4b8043e0/AirFoilCoordinatesArraysActivity.jl))*

### NACA 2412 Polar Plots Under Various Conditions

**Note:** As smoothing a polar plot requires extensive user input and trial and error, only the first plot/case will be extended and smoothed, the rest will simply be given as is.

**M = 0, Re = 10,000, Ncrit = 7**

<img title="" src="file:///C:/Users/child/Documents/Flow_Lab/Onboarding/XfoilProject/OutputImages/NACA2412PlotRe10000M0N7.jpg" alt="NACA2412Polar" data-align="center" width="364">

*Figure 2: Plot under the given conditions (M = 0, Re = 10,000, Ncrit = 7) at different angles of attack given in degrees (Originally -1 to 10 deg before extending). All the plots were smoothed, and the Cm and Cl plots were extended and smoothed from -pi to pi. For methods, see the source code linked in the appendix.*

**M = 0, Re = 100000, Ncrit = 7**

<img title="" src="file:///C:/Users/child/Documents/Flow_Lab/Onboarding/XfoilProject/OutputImages/NACA2412PlotRe100000M0N7.jpg" alt="NACA2412Polar" data-align="center" width="323">

*Figure 3: Plot under the given conditions (M = 0, Re = 100000, Ncrit = 7) at different angles of attack given in degrees (from -10 to 15 deg)*

**M = 0, Re = 1000000, Ncrit = 7**

<img title="" src="file:///C:/Users/child/Documents/Flow_Lab/Onboarding/XfoilProject/OutputImages/NACA2412PlotRe1000000M0N7.jpg" alt="FlowLab_Onboarding/NACA2412PlotRe1000000M0N7.jpg at f6e43fcf5526f132a274fbd0b73e05784ecc6aa7 · JacobChild/FlowLab_Onboarding · GitHub" data-align="center" width="334">

*Figure 4: Plot under the given conditions (M = 0, Re = 1000000, Ncrit = 7) at different angles of attack given in degrees (from -10 to 15 deg)*

**M = 0.02, Re = 10000, Ncrit = 11**

<img title="" src="file:///C:/Users/child/Documents/Flow_Lab/Onboarding/XfoilProject/OutputImages/NACA2412PlotRe10000M02N11.jpg" alt="FlowLab_Onboarding/NACA2412PlotRe10000M02N11.jpg at f6e43fcf5526f132a274fbd0b73e05784ecc6aa7 · JacobChild/FlowLab_Onboarding · GitHub" data-align="center" width="393">

*Figure 5: Plot under the given conditions (M = 0.02, Re = 10000, Ncrit = 11) at different angles of attack given in degrees (from -10 to 15 deg)*

**M = 0.02, Re = 100000, Ncrit = 11**

<img title="" src="file:///C:/Users/child/Documents/Flow_Lab/Onboarding/XfoilProject/OutputImages/NACA2412PlotRe100000M02N11.jpg" alt="Plot" data-align="center" width="383">

*Figure 6: Plot under the given conditions (M = 0.02, Re = 100000, Ncrit = 11) at different angles of attack given in degrees (from -10 to 15 deg). It should be noted that my code failed under these conditions when not using the same coordinates as Xfoil, ie, it was necessary to have more coordinate points than the default settings in PolarPlotter.jl to have converged results.*

**M = 0.02, Re = 1000000, Ncrit = 11**

<img title="" src="file:///C:/Users/child/Documents/Flow_Lab/Onboarding/XfoilProject/OutputImages/NACA2412PlotRe1000000M02N11.jpg" alt="FlowLab_Onboarding/NACA2412PlotRe1000000M02N11.jpg at f6e43fcf5526f132a274fbd0b73e05784ecc6aa7 · JacobChild/FlowLab_Onboarding · GitHub" data-align="center" width="371">

*Figure 7: Plot under the given conditions (M = 0.02, Re = 1000000, Ncrit = 11) at different angles of attack given in degrees (from -10 to 15 deg)*

### 

### Validation and Verification through Comparison of Results

*This was done between Xfoil, my adapted Xfoil.jl code, Airfoiltools.com, and experimental data.*

**Validation with known theoretical data**

- Example From Airfoiltools.com 
  
  - M = 0, Ncrit = 9, Re = 50,000 , iter = 500 (their iteration limit is unknown)
  
  - Using downloaded coordinate file from airfoiltools.com- NACA2412.dat
  
  - *Note: Xfoil.jl uses 140 panels by default, Xfoil uses 160, when I updated the panels it matched*
  
  *Table 1: Coefficient of lift comparison between Xfoil, My code, and Airfoiltools.com at various angles of attack*

| Angle | Xfoil  | MyCode  | Online |
|:-----:|:------:|:-------:|:------:|
| 0     | -.0589 | -.05884 | -.0589 |
| 1     | .1375  | .13749  | .1375  |
| 10    | 1.1591 | 1.15908 | 1.1591 |

**Error Analysis**

```matlab
%Error calculations at the 0 degree angle of attack
% at 0 deg the calculations appear to be the furthest off, so this
%will give us the worst case scenario
m = -.05884;     %the value from my code at 0 degs
x = -.0589;   %the value from xfoil at 0 degs

RelativeError = (m-x)/x*100 => -0.10186757%
```

**Validation Error Discussion**

Using the same coordinate file, and changing the number of panes in my adapted xfoil.jl code led to very little relative error between the Xfoil output and my adapted Xfoil.jl output. Double checking the Xfoil results with online validated that the test was being run with all of the correct inputs. While there is very little error, it is important to discuss the source of what error there is. My code (ie the Julia wrapper) uses the exact same solver as Xfoil but has some additional helps built in that could cause differences. After discussion with Taylor McDonald I learned that Xfoil.jl has additional methods built in to help the solver converge to a solution (ie "percussive maintenance" etc). Further inspection of my code and results revealed that it converged on all given output points, however Xfoil did not. Thus the source of the "error" was revealed and can be safely assumed to be within allowed tolerances with a relative error of less than a percent.

**Verification through comparison with experimental data**

The experimental data and plots that were used can be found [here](https://arc.aiaa.org/doi/pdfplus/10.2514/6.2018-1277.c1)

- As the experimental data was not tabulated, but was instead plotted, the estimated stall angle and Coefficient of lift at that angle will be reported and compared

*Table 2: Comparison of the estimated stall angles and coefficients of lift under these conditions: Re = 3,100,000, M = 0, and Ncrit = 9. (Note a mach number was not specified, so M = 0 was assumed and used)*

|                                              | XFoil  | MyCode | Experimental |
| -------------------------------------------- | ------ | ------ | ------------ |
| Estimated Stall Angle (deg)                  | 18     | 18     | 15           |
| Coefficient of Lift at estimated stall angle | 1.7689 | 1.7733 | 1.6          |

| ![myplot](C:\Users\child\Documents\Flow_Lab\Onboarding\XfoilProject\OutputImages\NACA2412PlotRe3100000M0N9.jpg) | ![expplot](C:\Users\child\Documents\Flow_Lab\Onboarding\XfoilProject\OutputImages\ExperimentalNACA2412Re3100000M0N9.jpg) |
|:---------------------------------------------------------------------------------------------------------------:|:------------------------------------------------------------------------------------------------------------------------:|
| *Figure 8: Theoretical output calculated from the PolarPlotter.jl code*                                         | *Figure 9: Experimental output found in a wind tunnel and cited in the link above*                                       |

**Error Analysis**

```matlab
%Relative error analysis 
m = 18;     %predicted from "my code" and xfoil
e = 15;    %the experimental stall angle as shown in the plot

%Stall angle error
RelativeStallAngleError = (m - e) / e * 100
= 20.0%

%Cl error- m and e were redefined with the appropriate Cl values
RelativeClError = (m - e) / e * 100
= 10.83%
```

**Verification Error Discussion**

The error between the theoretical and experimental data is quite large. It must be remembered however, that the experimental data points were estimated off of a plot, it was assumed that M = 0, and there were unknown assumptions made by the  original researchers. It follows that the aggregation of these error sources would lead to a quite large error overall. The general trends that we see with the theoretical data match what we would expect to see- the theoretical model performs better (not stalling until later), than the real life model, and the linear region and shape of the curve also generally match between cases. The theoretical model overpredicting performance of the airfoil can be quite dangerous as expecting the airfoil to perform one way and it falling short could lead to problems. 

### Overall Discussion

**Reynolds Number Trends**

The Mach number and Ncrit values were left constant as the Reynolds number was varied under two different sets of conditions. A high Reynolds number means that the effects of inertia are greater then the effects of the viscosity of the fluid. In simple terms this means the flow / boundary layer will become turbulent quicker. In the case of airfoils, turbulent boundary layers keep the flow attached to the airfoil, thus delaying flow separation longer at higher Reynolds numbers.  When comparing the coefficients of lift seen under the first set of conditions (for example), this effect can be seen in the higher Cl values. While turbulence and delayed flow separation play a part in increased lift, higher flow velocity is a larger cause of the increase in lift. With the fluid properties kept the same, a higher Reynolds number comes from a higher fluid velocity. For an angle of attack of 10 degrees, Cl at Re = 10,000 is apx .4 and at Re = 100000 it is apx 1.2, and at Re = 1000000 it is apx 1.5. This shows that higher Reynolds numbers allow for the generation of more lift due to the flow velocity increasing and that same flow staying attached to the airfoil longer. As the Reynolds number increases the slope (or rate of that increase) of the coefficient of lift also becomes greater, meaning that the coefficient of lift increases more per degree at higher Reynolds numbers than lower numbers.

These effects are even more apparent when looking at the Coefficient of Drag. A large portion of drag comes from "pressure drag", or the drag that comes because of the pressure differential between the front and trailing edge of the airfoil. This pressure differential is caused by flow separation, so it follows that if there is less flow separation there is less drag. These results are seen when comparing the Cd plots under the second condition at a 10 degree angle of attack (for example). At Re = 10000, Cd is approximately 1.1, at Re = 100000 it is approximately .04, and at Re = 1000000 it is approximately .016. Thus, it can be seen that drag is drastically reduced with an increase of the Reynolds number.

**Angle of Attack Trends**

As the angle of attack was varied from negative (nose down) to positive (nose up), regardless of the Reynolds number and other conditions, a few general trends could be seen. The Coefficient of lift generally grows more positive as the angle of attack is increased until it eventually reaches stall (where the airfoil stops producing as much lift). It should be noted that at some negative angles of attack the Coefficient of lift briefly stops growing more positive (at this point Cl is already negative) but has a negative slope for a few degrees before reversing again and becoming more positive as seen in Figure 10. This is due to the streamline/flow being largely unattached from the airfoil, and then reattaching as the negative angle of attack is less extreme and briefly generating more negative lift. 

<img title="" src="file:///C:/Users/child/Documents/Flow_Lab/Onboarding/XfoilProject/OutputImages/NegativeClExample.jpg" alt="NegativeCl" width="246" data-align="center">

*Figure 10: Coefficient of lift with the slopes outlined at negative angles of attack. Taken from the M = 0.02, Re = 1000000, Ncrit = 11 Condition*

The Coefficient of Drag displays more intuitive behavior in that as the angle of attack travels from negative to positive angle at which the least drag occurs is typically within a few degrees of 0. This makes sense as that is when the airfoil would be the most streamlined and the flow would be the least separated, thus generating the least friction and pressure drag respectively. Keep in mind that as was already shown above the effects of pressure drag are much larger than friction drag, and thus the angles of attack at which flow separation is larger are where the Coefficient of Drag is also much larger.

The Coefficient of Moment in all cases stays negative within our range of angles of attack. A negative Coefficient of Moment means that the airfoil is wanting to pitch nose down. This is desirable behavior as it means that the airfoil will resist stall (the angle at which an airfoil no longer produces more lift). The moment Coefficient is closely related to the other coefficients. In general both Lift and Drag give a positive moment to the airfoil (ie cause the nose to pitch up). This can be explained by the pressure distribution along the airfoil. At a higher angle of attack, both drag and lift will be trying to rotate the airfoil nose up, but a closer look at the pressure distribution shows that the locations and magnitudes of the pressure actually cause a negative moment coefficient. At higher Reynolds numbers as the angle of attack increases, the coefficient of moment generally becomes less negative, and at lower Reynolds numbers as the angle of attack increases, the coefficient of moment generally becomes more negative.

<img title="" src="file:///C:/Users/child/Documents/Flow_Lab/Onboarding/XfoilProject/OutputImages/StallExample.jpg" alt="StallExample" data-align="center" width="390">

*Figure 11: Example of the location of positive stall as outlined in the different plots. Taken from the M = 0.02, Re = 1000000, Ncrit = 11 Condition*

**Stall Trends**

When stall occurs the airfoil generates less and less lift rather than more than before (for positive stall). This typically happens at higher angles off attack, with relatively high (or rapidly increasing) coefficient of drag, and a very small moment coefficient (ie not a strong desire to pitch the nose down). These trends can be seen in figure 11. Stall can also occur in the negative Coefficient of Lift region. This can happen when the airfoil is beginning to generate more positive lift (ie still negative lift but with a positive slope) and then suddenly becoming more negative again as previously shown in figure 10. Negative stall also occurs when Cl is negative, but typically comes from the airfoil being inverted (extremely high or low angles of attack), and then just rotating out of its inversion as shown in figure 12.

<img title="" src="file:///C:/Users/child/Documents/Flow_Lab/Onboarding/XfoilProject/OutputImages/NegativeStallExample.jpg" alt="NegativeStallExample" data-align="center" width="290">

*Figure 12: Negative stall example as shown on the extended plot. Taken from the M = 0.0, Re = 10000, Ncrit = 7 Condition*

### Conclusion

A much deeper understanding of the Coefficients of Lift, Drag, and Moment were gained, as well as a more intuitive grasp of how they vary with angle of attack. A proficiency at basic airfoil analysis, plotting, and comparison was gained, but with XFoil, and Xfoil.jl. Julia coding skills were improved with an increased knowledge of how to use functions in different files and how to use different math packages. The code was both validated and verified and will be safe to use (within normal conditions) in the future.

#### Appendix

**A. High Level Coding Approach**

- Figure 1 was generated through the AirFoilCoordinatesFunction.jl file, with the lines to make a function commented out and the lines to plot uncommented. The same output can come from the AirFoilCoordinatesArraysActivity.jl file with no changes as that was made to stand alone. Both plots are validated.

- Figure 2 was generated through the ExtenderNSmoother.jl file, which called on the PolarPlotter.jl function to generate the coefficient values to be extended and smoothed. PolarPlotter.jl called on the AirFoilCoordinatesFunction.jl function to generate the airfoil coordinates to be fed to the xfoil solver in PolarPlotter. The resulting data from PolarPlotter.jl was smoothed by creating an overlayed 1-D spline and evaluating it. That was made possible through the Dierckx package. The smoothed Cl and Cd data was then extended through the Viterna function from -pi to pi, made possible by the CCBlade package. The extended data was then smoothed and plotted.

- Figures 3-7 (ie non extended and smoothed plots) were generated through the TempPolarPlotter.jl file which was virtually identical to PolarPlotter.jl, but is not a function and was made to output the needed plots and data slightly differently than ExtenderNSmoother.jl did.

**B. Source Code**

[ExtenderNSmoother.jl](https://github.com/JacobChild/FlowLab_Onboarding/blob/00d1b4217b55f39b4ca7fdd2090e0cdfc3032723/XfoilProject/ExtenderNSmoother.jl)

[PolarPlotter.jl](https://github.com/JacobChild/FlowLab_Onboarding/blob/00d1b4217b55f39b4ca7fdd2090e0cdfc3032723/XfoilProject/PolarPlotter.jl)

[AirFoilCoordinatesFunction.jl](https://github.com/JacobChild/FlowLab_Onboarding/blob/00d1b4217b55f39b4ca7fdd2090e0cdfc3032723/XfoilProject/AirFoilCoordinatesFunction.jl)

[TempPolarPlotter.jl](https://github.com/JacobChild/FlowLab_Onboarding/blob/00d1b4217b55f39b4ca7fdd2090e0cdfc3032723/XfoilProject/TempPolarPlotter.jl)

For all other code, output images, and this document (in Markdown) see the [XfoilProject](https://github.com/JacobChild/FlowLab_Onboarding/tree/XFoil-Project-Branch/XfoilProject) folder on my GitHub. 

*Note: All of these links currently take the user to a branch off of main. Once everything in the project has been validated and verified, a pull request will be created and this branch will be merged with Main. It is yet to be seen if the links will work once that happens.*
