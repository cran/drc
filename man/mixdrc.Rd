\name{mixdrc}

\alias{mixdrc}

\title{Fitting multiple dose-response curves with random effects}

\description{
  'mixdrc' fits non-linear mixed models where random effects are assigned to some or all of the parameters.
}

\usage{
mixdrc(object, random, data, startVal)
}

\arguments{
  \item{object}{an object of class 'drc'.}
  \item{random}{a character string specifying the random effects. Specification is as in \code{\link[nlme]{nlme}}.}
  \item{data}{a data frame containing the variables in the model (not optional).}
  \item{startVal}{a numeric vector of initial parameter values for the fixed effects (optional).}
}

\details{
  The function searches through a range of initial values in order to obtain convergence of the estimation procedure.
}

\value{
  An object of classes 'drc' and 'nlme'. 
}

% \references{ ~put references to the literature/web site here ~ }

\author{Christian Ritz}

\note{
  This functionality will be completely revamped and extended during 2009.
}

\seealso{Non-linear mixed models can be fitted using \code{\link[nlme]{nlme}}.}

\examples{

## First fitting a model without random effects
## Oops: ideally a Box-Cox transformation should be used
spinach.m1 <- drm(SLOPE~DOSE, CURVE, data = spinach, fct = LL.4(),
pmodels = data.frame(HERBICIDE, 1, 1, HERBICIDE))

## Then fitting the same model (the same fixed effects)
##  but with random effects (assigned to d)

spinach.m2 <- mixdrc(spinach.m1, random="d~1|CURVE", data = spinach)
compParm(spinach.m2, "e", "-")

## Fitting a model without random effects
##  but with a Box-Cox transformation

spinach.m3 <- drm(SLOPE~DOSE, CURVE, data = spinach, fct = LL.4(),
pmodels = data.frame(HERBICIDE, 1, 1, HERBICIDE))

## Then fitting the same model with random effects assigned to d

spinach.m4 <- mixdrc(spinach.m3, random="d~1|CURVE", data = spinach)
compParm(spinach.m4, "e", "-")

## Fitting a model without random effects, but fixing the c parameter at 0
spinach.m5 <- drm(SLOPE~DOSE, CURVE, data = spinach, fct = LL.3(),
pmodels = data.frame(HERBICIDE, 1, HERBICIDE))

## Fitting the corresponding mixed model with random effects on d
spinach.m6 <- mixdrc(spinach.m5, random="d~1|CURVE", data = spinach)

}
\keyword{models}
\keyword{nonlinear}
