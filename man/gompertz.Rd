\name{gompertz}

\alias{gompertz}

\title{The Gompertz (Weibull) function}

\description{
  'gompertz' provides a very general way of specifying Gompertz/Weibull dose response functions,
  under various constraints on the parameters.
}

\usage{
gompertz(lowerc = c(-Inf, -Inf, -Inf, -Inf),
         upperc = c(Inf, Inf, Inf, Inf), 
         fixed = c(NA, NA, NA, NA), 
         names = c("b", "c", "d", "e"), useDer = FALSE)
}

\arguments{
  \item{lowerc}{numeric vector. The lower bound on parameters. Default is minus infinity.}
  \item{upperc}{numeric vector. The upper bound on parameters. Default is plus infinity.}
  \item{fixed}{numeric vector. Specifies which parameters are fixed and at what value they are fixed. NAs for parameter that are not fixed.}
  \item{names}{a vector of character strings giving the names of the parameters (should not contain ":"). The default is reasonable (see under 'Usage').
               The order of the parameters is: b, c, d, e (see under 'Details').}
  \item{useDer}{logical. If TRUE derivatives are supplied, otherwise they are not supplied. Not yet implemented!}
}

\details{
  The default arguments yield a four-parameter Gompertz model
  \deqn{ f(x) = c + (d-c) \exp(-\exp(b(\log(x)-e))).}
  
  The function is asymmetric about the inflection point, that is the parameter \eqn{\exp(e)}.
  
  Due to the parameterisation used, fitting a Gompertz model may be more robust
  than fitting a logistic model.
}

\value{
  The value returned by the 'gompertz' is a list with the following components
  \item{fct}{The dose response function.}
  \item{ssfct}{The self starter function.}
  \item{deriv1}{The first derivative.}  
  \item{deriv2}{The second derivative.}  
  \item{lowerc}{The lower bounds on the parameters.}
  \item{upperc}{The upper bounds on the parameters.}
  \item{edfct}{The ED function.}
  \item{sifct}{The SI function.}    
}

\references{Seber, G. A. F. and Wild, C. J (1989) \emph{Nonlinear Regression}, New York: Wiley \& Sons (pp. 330--331).}

\author{Christian Ritz}

\note{This function is for use with the function \code{\link{multdrc}}.}

\seealso{Special cases of the function 'gompertz' are \code{\link{g3}} and \code{\link{g4}}.}

\examples{

## A four-parameter Gompertz model with b fixed at 1
model1 <- multdrc(FA, fct=gompertz(fixed=c(1, NA, NA, NA)))
summary(model1)

## A four-parameter Gompertz model with the constraint b>3
model2 <- multdrc(FA, fct=gompertz(lowerc=c(3, -Inf, -Inf, -Inf)), control=mdControl(constr=TRUE))
summary(model2)

rm(model1, model2)

}
\keyword{models}
\keyword{nonlinear}
