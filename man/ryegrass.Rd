\name{ryegrass}

\alias{ryegrass}

\docType{data}

\title{Effect of ferulic acid on growth of ryegrass}

\description{
  A single dose-response curve.
}

\usage{data(ryegrass)}

\format{
  A data frame with 24 observations on the following 2 variables.
  \describe{
    \item{rootl}{a numeric vector (response)}  
    \item{conc}{a numeric vector (predictor)}
  }
}

\details{
  The data are part of a study to investigate the joint action
  of phenolic acids on root growth inhibition of perennial ryegrass (\emph{Lolium perenne L}).

  \code{conc} is the concentration of ferulic acid is in mM, and \code{rootl} is the root length 
  of perennial ryegrass. 
}

\source{
  Inderjit and J. C. Streibig, and M. Olofsdotter (2002) Joint action of
  phenolic acid mixtures and its significance in allelopathy
  research, \emph{Physiologia Plantarum}, \bold{114}, 422--428, 2002.
}

\examples{

## Displaying the data set
ryegrass

## Fitting a four-parameter Weibull model (type 2)
ryegrass.m1 <- drm(rootl ~ conc, data = ryegrass, fct = W2.4())

## Displaying a summary of the model fit
summary(ryegrass.m1)

## Plotting the fitted curve together with the original data
plot(ryegrass.m1)

## Fitting a four-parameter Weibull model (type 1)
ryegrass.m2 <- drm(rootl ~ conc, data = ryegrass, fct = W1.4())
plot(ryegrass.m2)

## Fitting a four-parameter log-logistic model
## with user-defined parameter names
ryegrass.m3 <- drm(rootl ~ conc, data = ryegrass, 
fct = LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))
summary(ryegrass.m3)

}
\keyword{datasets}

