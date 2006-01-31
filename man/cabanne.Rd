\name{cabanne}
\alias{cabanne}

\docType{data}

\title{data name/kind}

\description{
  A concise (1-5 lines) description of the dataset.
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
  reference to a publication or URL from which the data were obtained
}

\references{
  possibly secondary sources and usages
}

\examples{

model1 <- multdrc(drymat~dose, curve, data=cabanne, cm="0")
plot(fitted(model1), residuals(model1))

## Fitting one curve with variance as a power of the mean
model2 <- multdrc(drymat~dose, data=cabanne[c(1:40,141:240),], varPower=TRUE)

rm(model1, model2)

}
\keyword{datasets}
