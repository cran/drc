\name{mselect}

\alias{mselect}

\title{Model selection}

\description{
  Model selection by comparison of different models using the maximum log likelihood value, 
  Akaike's information criterion (AIC), the estimated residual variance 
  and the p-value from a lack-of-fit test as criteria.
}

\usage{
  mselect(object, fctList = NULL, nested = FALSE)
}

\arguments{
  \item{object}{an object of class 'drc'.}
  \item{fctList}{a list of non-linear functions to be compared.}
  \item{nested}{logical. TRUE results in F tests between adjacent models (in 'fctList'). 
  Only sensible for nested models.}
}

\details{
  The maximum likelihood cannot be used for comparison unless the models are nested.

  For Akaike's information criterion and the residual variance: the smaller the better and 
  for lack-of-fit test: the larger the better. Note that the residual variance is only available for 
  continuous dose-response data.
}

\value{
  A matrix with one row for each model and one column for each criterion.
}

%\references{}

\author{Christian Ritz}

%\note{}

%\seealso{}

\examples{

### Example with continuous data
## Fitting initial four-parameter log-logistic model
ryegrass.m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

## Model selection
mselect(ryegrass.m1, list(LL.3(), LL.5(), W1.3(), W1.4(), W2.4(), baro5()))

mselect(ryegrass.m1, list(LL.5()), nested = TRUE)  # nested models

### Example with quantal data
## Fitting initial two-parameter log-logistic model
earthworms.m1 <- drm(number/total~dose, weights=total, 
data = earthworms, fct = LL.2(), type = "binomial")

## Comparing 4 models
mselect(earthworms.m1, list(W1.2(), W2.2(), LL.3()))
# model selection based AIC or p value not working in this example!!!

}
\keyword{models}
\keyword{nonlinear}
