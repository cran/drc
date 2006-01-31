\name{FA}
\alias{FA}
\docType{data}
\title{Effect of Ferulic Acid on growth of ryegrass}
\description{
  A single dose-response curve.
}
\usage{data(FA)}
\format{
  A data frame with 24 observations on the following 2 variables.
  \describe{
    \item{MEANLR}{a numeric vector (response)}  
    \item{MM}{a numeric vector (dose)}
  }
}
\details{
  The data are part of a study to investigate the joint action
  of phenolic acids on root growth inhibition of perennial ryegrass (\emph{Lolium perenne L}).

  MM is the concentration of ferulic acid in mM, and MEANLR is the root length of perennial ryegrass.} 
}
\source{

Inderjit and J. C. Streibig, and M. Olofsdotter (2002) Joint action of
phenolic acid mixtures and its significance in allelopathy
research, \emph{Physiologia Plantarum}, \bold{114}, 422--428, 2002.
}
\examples{

## Displaying the data set
FA

## Fitting a four-parameter logistic regression 
## with user-defined parameter names
model1 <- multdrc(MEANLR~MM, data=FA, 
fct = l4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))

## Getting a summary of the fit
summary(model1)


## Fitting a model with very small dose values
FA2 <- FA
FA2[,2] <- FA2[,2]/1e6
model2 <- multdrc(MEANLR~MM, data=FA2)

## Removing the model objects
rm(model1, model2, FA2)

}
\keyword{datasets}
