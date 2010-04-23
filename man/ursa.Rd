\name{ursa}

\Rdversion{1.1}

\alias{ursa}

\title{
  Model function for the universal response surface approach (URSA) for the quantitative assessment of drug interaction
}

\description{
  URSA provides a parametric approach for modelling the joint action of several agents. The model allows quantification of synergistic effects through a single parameter.
}

\usage{
ursa(fixed = c(NA, NA, NA, NA, NA, NA, NA), names = c("b1", "b2", "c", "d", "e1", "e1", "f"), ssfct = NULL)
}

\arguments{
  \item{fixed}{numeric vector. Specifies which parameters are fixed and at what value they are fixed. NAs for parameter that are not fixed.}
  \item{names}{a vector of character strings giving the names of the parameters. The default is reasonable.}
  \item{ssfct}{a self starter function to be used (optional).}  
}

\details{
  The model function is defined implicitly through an appropriate equation. The details are found in Syracuse and Greco (1986), Greco et al (1990, 1995).
}

\value{
  A list containing the nonlinear function, the self starter function, and the parameter names.
}

\references{
  Finney, D. J. (1979) Bioassay and the Practise of Statistical Inference, 
  \emph{Int. Statist. Rev.}, \bold{47}, 1--12.

  Syracuse, K. C. and Greco, W. R. (1986) Comparison between the method of Chou and Talalay and a new method for the assessment of the combined effects
  of drugs: A Monte-Carlo simulation study, \emph{Proceedings of the Biopharmaceutical Section of the American Statistsical Association}, 127--132.

  Greco, W. R. and Park H. S. and Rustum, Y. M. (1990) Application of a New Approach for the Quantitation of Drug Synergism
  to the Combination of cis-Diamminedichloroplatinum and 1-beta-D-Arabinofuranosylcytosine, \emph{Cancer Research}, \bold{50}, 5318--5327.
  
  Greco, W. R. Bravo, G. and Parsons, J. C. (1995) The Search for Synergy: A Critical Review from a Response Surface Perspective,
  \emph{Pharmacological Reviews}, \bold{47}, Issue 2, 331--385.

}

\author{
  Christian Ritz after an idea by Hugo Ceulemans
}

%\note{}

\seealso{
  Other models for fitting mixture data are the Hewlett and Voelund models \code{\link{mixture}}, ~~~
}

\examples{

## Data from Greco et al (1995)
d1 <- c(0, 0, 0, 0, 0, 0, 0, 0, 2, 5, 10, 20, 50, 2, 2, 2, 2, 2, 5, 5, 5, 5, 5,10, 10, 10, 10, 10, 20, 20, 20, 20, 20, 50, 50, 50, 50, 50)

d2 <- c(0, 0, 0, 0.2, 0.5, 1, 2, 5, 0, 0, 0, 0, 0, 0.2, 0.5, 1, 2, 5, 0.2, 0.5, 1, 2, 5, 0.2, 0.5, 1, 2, 5, 0.2, 0.5, 1, 2, 5, 0.2, 0.5, 1, 2, 5)

effect <- c(106.00, 99.20, 115.00, 79.20, 70.10, 49.00, 21.00, 3.83, 74.20, 71.50, 48.10, 30.90, 16.30, 76.30, 48.80, 44.50, 15.50, 3.21, 56.70, 47.50,
26.80, 16.90, 3.25, 46.70, 35.60, 21.50, 11.10, 2.94, 24.80, 21.60, 17.30, 7.78, 1.84, 13.60, 11.10, 6.43, 3.34, 0.89)

greco <- data.frame(d1, d2, effect)

## Fitting the URSA model
greco.m1<-drm(effect~d1, data=greco, fct=ursa(fixed=c(NA,NA,0,NA,NA,NA,NA)), pmodels=list(~1,~1,~1,~I(1/d1)-1,~I(1/d2)-1,~1))

plot(fitted(greco.m1), residuals(greco.m1))
summary(greco.m1)

## Fitting the URSA model using values from Greco et al (1995) (p. 364)
greco.m2<-drm(effect~d1, data=greco, fct=ursa(fixed=c(NA,NA,0,NA,NA,NA,NA)), pmodels=list(~1,~1,~1,~I(1/d1)-1,~I(1/d2)-1,~1), start=c(-1.05, -2.04,95.1,11.1,1.07,0.52))
# same fit as above

## Adjusting for variance heterogeneity
greco.m3<-boxcox(greco.m2, method="anova")

plot(fitted(greco.m3), residuals(greco.m3))  # improved, not great
summary(greco.m3)

## Getting closer to Greco's estimates
#greco.m4<-drm(effect~d1, data=greco, fct=ursa(fixed=c(NA,NA,0,NA,NA,NA,NA)), pmodels=list(~1,~1,~1,~I(1/d1)-1,~I(1/d2)-1,~1), bcVal=0)  # not perfect
greco.m4<-drm(effect~d1, data=greco, fct=ursa(fixed=c(NA,NA,0,NA,NA,NA,NA)), pmodels=list(~1,~1,~1,~I(1/d1)-1,~I(1/d2)-1,~1), bcVal=0, start=coef(greco.m3))

plot(fitted(greco.m4), residuals(greco.m4))  # better!

summary(greco.m4)
}

\keyword{models}
\keyword{nonlinear}

%\concept{Hill 4-parameter four-parameter}
