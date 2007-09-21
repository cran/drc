\name{heartrate}

\alias{heartrate}

\docType{data}

\title{Heart rate baroreflexes for rabbits}

\description{
  The dataset contains measurements of mean arterial pressure (mmHG) and heart rate (b/min) for a baroreflex curve.
}

\usage{data(heartrate)}

\format{
  A data frame with 18 observations on the following 2 variables.
  \describe{
    \item{\code{pressure}}{a numeric vector containing measurements of arterial pressure.}
    \item{\code{rate}}{a numeric vector containing measurements of heart rate.}
  }
}

\details{
  The dataset is an example of an asymmetric dose-response curve, that is not
  easily handling using the log-logistic or Weibull models, whereas the baro5 model
  provides a nice fit.
}

\source{Ricketts, J. H. and Head, G. A. (1999) A five-parameter logistic equation for investigating asymmetry of
curvature in baroreflex studies. \emph{Am. J. Physiol. (Regulatory
Integrative Comp. Physiol. 46)}, \bold{277}, 441--454.}

\examples{

## Fitting the baro5 model
model1 <- multdrc(rate~pressure, data=heartrate, fct=baro5())
plot(model1)

rm(model1)

}
\keyword{datasets}
