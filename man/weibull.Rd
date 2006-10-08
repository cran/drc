\name{weibull}

\alias{weibull}
\alias{weibull2}

\title{The Weibull function}

\description{
  'weibull' and 'weibull2' provide a very general way of specifying Weibull dose response functions,
  under various constraints on the parameters.
}

\usage{
  weibull(lowerc = c(-Inf, -Inf, -Inf, -Inf),
          upperc = c(Inf, Inf, Inf, Inf), 
          fixed = c(NA, NA, NA, NA), 
          names = c("b", "c", "d", "e"),
          scaleDose = TRUE, useDer = FALSE)
 
  weibull2(lowerc = c(-Inf, -Inf, -Inf, -Inf),
           upperc = c(Inf, Inf, Inf, Inf), 
           fixed = c(NA, NA, NA, NA), 
           names = c("b", "c", "d", "e"),
           scaleDose = TRUE, useDer = FALSE)
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
  \item{useDer}{logical. If TRUE derivatives are supplied, otherwise they are not supplied. Not yet implemented!}
}

\details{
  As pointed out in Seber and Wild (1989), there exist two different parameterisations of the Weibull model. They
  do not yield the same fitted curve for a given dataset (see under Examples).

  One four-parameter Weibull model ('weibull') is
  \deqn{ f(x) = c + (d-c) \exp(-\exp(b(\log(x)-\log(e)))).}
 
  Another four-parameter Weibull model ('weibull2') is
  \deqn{ f(x) = c + (d-c) (1 - \exp(-\exp(b(\log(x)-\log(e))))).}  
  
  Both four-parameter functions are asymmetric about the inflection point, that is the parameter \eqn{\exp(e)}. 
}

\value{
  The value returned by the 'weibull' is a list with the following components
  \item{fct}{The dose response function.}
  \item{ssfct}{The self starter function.}
  \item{deriv1}{The first derivative.}  
  \item{deriv2}{The second derivative.}  
  \item{lowerc}{The lower bounds on the parameters.}
  \item{upperc}{The upper bounds on the parameters.}
  \item{edfct}{The function used to calculate ED values.}
  \item{sifct}{The function used to calculate SI values.}    
}

\references{
  Seber, G. A. F. and Wild, C. J (1989) \emph{Nonlinear Regression}, New York: Wiley \& Sons (pp. 338--339).
}

\author{Christian Ritz}

\note{The functions are for use with the function \code{\link{multdrc}}.}

\seealso{
  For convenience several special cases of the function 'weibull' are available: 
  \code{\link{w2}}, \code{\link{w3}} and \code{\link{w4}}.  
}

\examples{

## Fitting two different Weibull models
model1 <- multdrc(FA, fct = weibull())
plot(model1, conLevel=0.5, legend=FALSE)
model2 <- multdrc(FA, fct = weibull2())
plot(model2, conLevel=0.5, type="add", col=2)
# you could also look at the ED values to see the difference

## A four-parameter Weibull model with b fixed at 1
model3 <- multdrc(FA, fct = weibull(fixed=c(1, NA, NA, NA)))
summary(model3)

## A four-parameter Weibull model with the constraint b>3
model4 <- multdrc(FA, fct = weibull(lowerc=c(3, -Inf, -Inf, -Inf)), 
control=mdControl(constr=TRUE))
summary(model4)

rm(model1, model2, model3, model4)
}
\keyword{models}
\keyword{nonlinear}
