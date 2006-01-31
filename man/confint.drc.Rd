\name{confint.drc}

\alias{confint.drc}

\title{Confidence Intervals for model parameters}

\description{
  Computes confidence intervals for one or more parameters in a model of class 'drc'.
}

\usage{
confint.drc(object, parm, level = 0.95, ...)
}

\arguments{
  \item{object}{a model object of class 'drc'.}
  \item{parm}{a specification of which parameters are to be given
          confidence intervals, either a vector of numbers or a vector
          of names. If missing, all parameters are considered.}
  \item{level}{the confidence level required.}
  \item{\dots}{additional argument(s) for methods. Not used.}
}

\details{
  The confidence intervals are based on asymptotic normality.
}

\value{
  A matrix (or vector) with columns giving lower and upper confidence limits for each parameter. These will be labelled as
  (1-level)/2 and 1 - (1-level)/2 in % (by default 2.5% and 97.5%).
}

%\references{ ~put references to the literature/web site here ~ }

\author{Christian Ritz}

%\note{}

%\seealso{}

\examples{

## Confidence interval for a single parameter
model1 <- multdrc(FA)
confint(model1, "e")

## Confidence intervals for all parameters
confint(model1)

rm(model1)

}
\keyword{models}
\keyword{nonlinear}
