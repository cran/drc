\name{bread.drc}

\alias{bread.drc}
\alias{estfun.drc}

\title{Bread and meat for the sandwich}

\description{
  Bread and meat for the sandwich estimator of the variance-covariance.
}

\usage{
  bread.drc(x, ...)
  
  estfun.drc(x, ...)
}

\arguments{
  \item{x}{object of class \code{drc}}
  \item{\dots}{additional arguments. At the moment none are supported.}
}

\details{
  More explanations are found in Zeileis (2006).
}

\value{
  The unscaled hessian is returned by \code{bread.drc}, whereas \code{estfun.drc} 
  returns the estimating function evaluated at the data and the parameter estimates.
}

\references{
  Zeileis, A. (2006) Object-oriented Computation of Sandwich Estimators,
  \emph{J. Statist. Software}, \bold{16}, Issue 9.
}

\author{Christian Ritz}

%\note{}

%\seealso{For other applications see \code{\link[sandwich]{sandwich}}.}

\examples{

## The lines below requires that the packages
## 'lmtest' and 'sandwich' are installed
#library(lmtest)
#library(sandwich)

#ryegrass.m1<-drm(rootl~conc, data=ryegrass, fct=LL.4())

#coeftest(ryegrass.m1)
#coeftest(ryegrass.m1, vcov=sandwich)

}

\keyword{models}
\keyword{nonlinear}
