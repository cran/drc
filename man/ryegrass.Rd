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
  of perennial ryegrass.} 
}

\source{
  Inderjit and J. C. Streibig, and M. Olofsdotter (2002) Joint action of
  phenolic acid mixtures and its significance in allelopathy
  research, \emph{Physiologia Plantarum}, \bold{114}, 422--428, 2002.
}

\examples{

## Displaying the data set
ryegrass

## Fitting a four-parameter logistic regression 
## with user-defined parameter names
model1 <- drm(rootl ~ conc, data = ryegrass, 
fct = LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))

## Displaying a summary of the fit
summary(model1)

}
\keyword{datasets}

