\documentclass[12pt]{texmemo} % originally by Rob Oakes; adapted by Alice Chen and now Jacob Child

\memoto{ Adam Cardoza, Graduate Student in the BYU FLOW Lab}

\memofrom{Jacob Child}

\memore{Findings from the DJI II Trade Study}

\memodate{\today} % or \today

\memosection{497R (Andrew Ning)}

\begin{document}

\maketitle

\highlight(Intro of what I will be describing in the paper) ?

*there are too many plots, include a few and put the rest attached below?

A set of custom code heavily utilizing the CCBlade.jl package was created in order to analyze variations on the DJI II propeller (A drone propeller found on the Phantom model *verify* ) The propeller was imported, analyzed, and its geometric and flow properties  indecently varied to see the trade-offs and effects on performance. The results and findings are presented below, with a focus on a discussion of the design space rather than the best/optimized parameters. 

*When and why do I limit things and how do the parameters influence the design space!

RVar findings

CVar findings

PVar findings

In each of the above discuss the effects on thrust, power, and figure of merit

Extra discussion on the Surface plots

Attachments

*do I put extra plots and code etc here?

Section I Discussion of the results that came from varying the blade radius from *value* to *value*

Section II Discussion of the results that came from varying the blade chord

Section III Discussion of the results that came from varying the blade pitch

Inside of each section the following items are discussed:

\setlength{\parskip}{0pt} % no new line in list

\begin{itemize}

    \item Coefficient of Thrust trends

    \item Coefficient of power and Torque trends 

    \item Figure of Merit Trends 

    \item  Overall Discussion

\end{itemize}

\setlength{\parskip}{0.5\baselineskip plus 2pt} % reset

*#?what is the above line for and do?*



Finally, Section IV concludes this memo with the a very general discussion of the design space and recommendations and limitations of the varied parameters. 

*Note:* Following the 4 main sections is an appendix containing:

\setlength{\parskip}{0pt} % no new line in list

\begin{itemize}

    \item Discussion of Surface Plots and their meaning

    \item Github overview and high level code explanation

    \item Extra plots and figures (undiscussed in favor of brevity)

    \end{itemize}



\section{Findings from varying the blade radius}

The purpose of this section is share the findings from varying the blade radius and how it affects aerodynamic coefficients and their implications on blade performance. 

\textit{Coefficient of Thrust} 

*input figure here for the ie RVarCoefPlot* *note:* there is currently an error in this plot? the 5000 rpm plot works just fine

The coefficient of thrust is the amount of thrust the rotor produces normalized with respect to the rpm. The plot above and general analysis was done at a fixed rpm of 5000 (for varied rpm results see the appendix). The blade radius was varied from 5cm to 25cm. Understandably CT is very small at the very small beginning blade radii. It increases quite quickly reaching 75% of the CT at 25cm when the blade radius is still only about 13cm. This makes sense as this blade was designed for a smaller drone and its original design radius is 12cm. After the 12-13cm mark CT increases much more slowly. This shows that while increasing the blade radius does always increase lift it does so at a slower and slower rate. Thinking of the geometry of the rotor this makes sense. When just the radius is varied the thickness and chord quickly become comparatively very small. This means that while the lifting surface area does increase, it is still a very small area over all and as the blade radius increases the blade becomes more and more "skinny". 

\textit{Coefficients of Power and Torque} *refer them to the same figure above as that has them also on it.*

The CQ and CP can be discussed together and interchangeably as CP is simply $2\pi*CQ$ this means that it is just scaled. The shape of the CQ increase closely follows the CT plot, generally, as CT increases, so also does CQ. This is likely because CT is increasing because the radius is, as the radius increases the blade becomes more and more difficult to turn (the center of mass moves farther from the hub and the blade mass increases). An increased torque means that the power required to turn the motor will also increase. 

\textit{Opposition Procedures.} The USPTO currently offers several ways of

\section{PVar Findings}

In this section I flesh out three major issues with software patents.

\textit{Slow Process.} It takes near

\textit{Low Quality.} Granted software

\textit{Costly Opposition.} If it is hard to ensure that granted patents are generally of high q

\section{CVar Findings}

In this section I propose three policies t

\textit{Clearer Guidelines.} tware patents and the tug-of-war between the SCOTUS and the CAFC, patent examiners and administrative law judges at the US

\textit{Stronger Collaboration.} \cite{lee} More and more students are majoring in computer science, more and more people are getting involved in open-sou

\textit{Smarter System.} Modernizing the IT system

\section{Surface Discussion}

filler bla bla baa

\section{Conclusion: Overall trends and implications, further research}

What the USPTO can do is limited.

\newpage

\printbibliography

\end{document}
