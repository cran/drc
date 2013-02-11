\name{beetGrowth}

\alias{beetGrowth}

\docType{data}

\title{Example of plant growth curve from Covarelli and Onofri (1998)}

\description{
  Crop weight (grams per square meter of dry matter) was recorded at several times
  after plant emergence
}

\usage{data(beetGrowth)}

\format{
  A data frame with 186 observations on the following 3 variables.
  \describe{
    \item{\code{DAE}}{a numeric vector}
    \item{\code{weightInf}}{a numeric vector}
    \item{\code{weightFree}}{a numeric vector}
  }
}

%\details{}

\source{
   Covarelli G. and Onofri A., 1998. Effects of timing of weed removal and emergence in sugarbeet. 
   Proceedings 6th EWRS Mediterranea Symposium, Montpellier, 65-72.
}

%\references{}

\examples{

model <- drm(weightInf ~ DAE, data  = beetGrowth, fct=gompGrowth.1())
plot(model, log="")
summary(model)
}
\keyword{datasets}
