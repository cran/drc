\name{cabanne}
\alias{cabanne}

\docType{data}

\title{data name/kind}

\description{
  Small plants, growing in pots in a green house, were sprayed with techincal 
  grade phenmidipham herbicide either alone or in mixture with an ester of oleic acid. 
  The plants were allowed to grow in the green house 14 days
  after treatment; then the dry matter was measured per pot.
}

\usage{data(cabanne)}

\format{
  A data frame with 240 observations on the following 3 variables.
  \describe{
    \item{\code{dose}}{a numeric vector}
    \item{\code{drymat}}{a numeric vector}
    \item{\code{curve}}{a numeric vector}
  }
}

\details{
  If necessary, more details than the description above 
}

\source{
Cabanne, F, Gaudry,JC & Streibig, JC (1999) Influence of alkyl oleates on efficacy 
of phenmedipham applied as an acetone:water solution on Galium aparine. 
Weed Research 39: 57-67.

}

\references{  possibly secondary sources and usages
}

\examples{

model1 <- multdrc(drymat~dose, curve, data=cabanne)
plot(fitted(model1), residuals(model1))

## Fitting one curve with variance as a power of the mean
model2 <- multdrc(drymat~dose, data=cabanne[c(1:40,141:240),], varPower=TRUE)

rm(model1, model2)

}
\keyword{datasets}
