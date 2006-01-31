\name{daphnids}

\alias{daphnids}

\docType{data}

\title{Daphnia test}

\description{
  The number of immobile daphnids --in contrast to mobile daphnids-- out of a total of 20 daphnids was counted 
  for several concentrations of a toxic substance.
}

\usage{data(daphnids)}

\format{
  A data frame with 16 observations on the following 4 variables.
  \describe{
    \item{\code{dose}}{a numeric vector}
    \item{\code{no}}{a numeric vector}
    \item{\code{total}}{a numeric vector}
    \item{\code{time}}{a factor with levels \code{24h} \code{48h}}
  }
}

\details{
  The same daphnids were counted at 24h and later again at 48h.
}

\source{
  Nina Cedergreen, Royal Veterinary and Agricultural University, Denmark.
}

%\references{}

\examples{

## Fitting a model with different parameters
##  for different curves
model1.daphnids <- multdrc(no/total~dose, time, weights=total, 
data=daphnids, fct=l2(), type="binomial")


## Goodness-of-fit test
anova(model1.daphnids)


## Fitting a model with a common intercept parameter
model2.daphnids <- multdrc(no/total~dose, time, weights=total, 
data=daphnids, fct=l2(), type="binomial", 
collapse=list(~1, ~time))


rm(model1.daphnids, model2.daphnids)
}
\keyword{datasets}
