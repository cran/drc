\name{ecvam}

\alias{ecvam}

\docType{data}

\title{Continuous data from ECVAM}

\description{
  Continuous data provided by ECVAM from tests of compounds denoted "a", "b" and so on, within a specific test system. 
}

\usage{
  data(ecvam)
}

\format{
  Data frames with 54 observations and 26 times 54 observations, respectively, on the following 2 variables.
  \describe{
    \item{\code{conc}}{a numeric vector of concentrations.}
    \item{\code{response}}{a numeric vector of response values.}
    \item{\code{experi}}{a character vector referring to the experiment (and compound).}
  }
}

\details{
  All datasets stem from experiments carried out under standard conditions, having replicates for all concentrations.
}

\source{
  Data are kindly provided by European Centre for the Validation of Alternative Methods (ECVAM), 
  EU Joint Research Centre, Ispra, Italy.
}

%\references{}

%\examples{}

\keyword{datasets}
