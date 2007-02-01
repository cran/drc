\name{mselect}

\alias{mselect}

\title{Model selection}

\description{
  Model selection by comparison of different models using 1) the maximum log likelihood value, 
  2) Akaike's information criteria (AIC), 3) the estimated residual variance 
  and 4) the p-value from a lack-of-fit test as criteria.
}

\usage{
mselect(object, fctList = NULL, nested = FALSE)
}

\arguments{
  \item{object}{an object of class 'drc'.}
  \item{fctList}{a list of non-linear functions to be compared.}
  \item{nested}{logical. TRUE results in F tests between adjacent (in 'fctList') models. 
  Only sensible for nested models.}
}

\details{
  Criterion 1) cannot be used for comparison unless the models are nested ().

  For criteria 2) and 3): the smaller the better; for criterion 4) the larger the better.
}

\value{
  A matrix with one row for each model and one column for each criterion.
}

%\references{}

\author{Christian Ritz}

%\note{}

%\seealso{}

\examples{

## Fitting initial model
m1 <- multdrc(ryegrass)

## Model selection
mselect(m1, list(l3(), l5(), W1.3(), W1.4(), weibull2(), baro5(), richards()))

mselect(m1, list(l5()))  # nested models

}
\keyword{models}
\keyword{nonlinear}
