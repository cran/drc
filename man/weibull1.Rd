\name{weibull1}

\alias{weibull1}
\alias{weibull2}

\title{The Weibull functions}

\description{
  'weibull' and 'weibull2' provide a very general way of specifying Weibull dose response functions,
  under various constraints on the parameters.
}

\usage{
  weibull1(lowerc = c(-Inf, -Inf, -Inf, -Inf),
           upperc = c(Inf, Inf, Inf, Inf), 
           fixed = c(NA, NA, NA, NA), 
           names = c("b", "c", "d", "e"),
           scaleDose = TRUE, fctName, fctText)
 
  weibull2(lowerc = c(-Inf, -Inf, -Inf, -Inf),
           upperc = c(Inf, Inf, Inf, Inf), 
           fixed = c(NA, NA, NA, NA), 
           names = c("b", "c", "d", "e"),
           scaleDose = TRUE, fctName, fctText)
}

\arguments{
  \item{lowerc}{numeric vector. The lower bound on parameters. Default is minus infinity.}
  \item{upperc}{numeric vector. The upper bound on parameters. Default is plus infinity.}
  \item{fixed}{numeric vector. Specifies which parameters are fixed and at what value they are fixed. 
  NAs for parameter that are not fixed.}
  \item{names}{a vector of character strings giving the names of the parameters (should not contain ":"). 
  The default is reasonable (see under 'Usage'). The order of the parameters is: b, c, d, e (see under 'Details').}
  \item{scaleDose}{logical. If TRUE dose values are scaled around 1 during estimation; 
                   this is required for datasets where all dose values are small. } 
  \item{fctName}{optional character string used internally by convenience functions.}  
  \item{fctText}{optional character string used internally by convenience functions.}                      
}

\details{
  As pointed out in Seber and Wild (1989), there exist two different parameterisations of the Weibull model. They
  do not yield the same fitted curve for a given dataset (see under Examples).

  One four-parameter Weibull model ('weibull1') is
  \deqn{ f(x) = c + (d-c) \exp(-\exp(b(\log(x)-\log(e)))).}
 
  Another four-parameter Weibull model ('weibull2') is
  \deqn{ f(x) = c + (d-c) (1 - \exp(-\exp(b(\log(x)-\log(e))))).}  
  
  Both four-parameter functions are asymmetric with inflection point at the dose \eqn{e}. 
}

\value{
  The value returned is a list containing the non-linear function, the self starter function
  and the parameter names.  
}

\references{
  Seber, G. A. F. and Wild, C. J (1989) \emph{Nonlinear Regression}, New York: Wiley \& Sons (pp. 338--339).
}

\author{Christian Ritz}

\note{The functions are for use with the function \code{\link{multdrc}}.}

\seealso{
  For convenience several special cases of the function 'weibull1' are available: 
  \code{\link{W1.2}}, \code{\link{W1.3}} and \code{\link{W1.4}}.  

  Special cases of 'weibull2' are: 
  \code{\link{W2.2}}, \code{\link{W2.3}} and \code{\link{W2.4}}.  
}

\examples{

## Fitting two different Weibull models
ryegrass.m1 <- drm(ryegrass, fct = weibull1())
plot(ryegrass.m1, conLevel=0.5)

ryegrass.m2 <- drm(ryegrass, fct = weibull2())
plot(ryegrass.m2, conLevel=0.5, add = TRUE, type = "fit", col=2)
# you could also look at the ED values to see the difference

## A four-parameter Weibull model with b fixed at 1
ryegrass.m3 <- drm(ryegrass, fct = weibull1(fixed=c(1, NA, NA, NA)))
summary(ryegrass.m3)

## A four-parameter Weibull model with the constraint b>3
ryegrass.m4 <- drm(ryegrass, fct = weibull1(), lowerl = c(3, -Inf, -Inf, -Inf), 
control=drmc(constr=TRUE))
summary(ryegrass.m4)

}
\keyword{models}
\keyword{nonlinear}
