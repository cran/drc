\name{TM}
\alias{TM}

\docType{data}

\title{A suite of bioassays}

\description{
  The data set consists of observations from 7 bioassays plus additional control measurements (dose equal to 0).
}

\usage{data(TM)}

\format{
  A data frame with 180 observations on the following 3 variables:
  \describe{
    \item{dose}{a numeric vector.}
    \item{pct}{a numeric vector.}
    \item{rgr}{a numeric vector.}
  }
}

\details{
  dose is the concentration of herbicide.
  pct is the grouping variable, separating observations into curves.
  rgr is the response (relative growth rate).  
}

\source{
  Cedergreen N. (2004) Unpublished bioassay data.
}

\note{A four-parameter Gompertz model provides an acceptable fit to the data.}

\examples{

TM[1:20,]  # displaying first 20 lines of the data set

## Fitting a Weibull model
## cControl measurements are specified using the 'cm' argument 
model1 <- multdrc(rgr~dose, pct, cm="999", data=TM, fct=g4())
summary(model1)

rm(model1)

}
\keyword{datasets}
