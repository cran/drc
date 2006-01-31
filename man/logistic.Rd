\name{logistic}

\alias{logistic}

\title{The logistic function}

\description{
  'logistic' provides a very general way of specifying logistic dose response functions,
  under various constraints on the parameters.
}

\usage{
  logistic(lowerc = c(-Inf, -Inf, -Inf, -Inf, -Inf), 
  upperc = c(Inf, Inf, Inf, Inf, Inf), fixed = c(NA, NA, NA, NA, NA), 
  names = c("b", "c", "d", "e", "f"), scaleDose = TRUE, useDer = FALSE)
}

\arguments{
  \item{lowerc}{numeric vector. The lower bound on parameters. Default is minus infinity.}
  \item{upperc}{numeric vector. The upper bound on parameters. Default is plus infinity.}
  \item{fixed}{numeric vector. Specifies which parameters are fixed and at what value they are fixed. 
               NAs for parameter that are not fixed.}
  \item{names}{a vector of character strings giving the names of the parameters (should not contain ":").
               The default is reasonable (see under 'Usage'). The order of the parameters is: b, c, d, e, f 
               (see under 'Details').}
  \item{scaleDose}{logical. If TRUE dose values are scaled around 1 during estimation; 
                   this is required for datasets where all dose values are small. }
  \item{useDer}{logical. If TRUE derivatives are supplied, otherwise they are not supplied. Not yet implemented!}
}

\details{
  The default arguments yields the five-parameter logistic function given by the expression
  \deqn{ f(x) = c + \frac{d-c}{(1+\exp(b(\log(x)-\log(e))))^f}}
  
  If the parameter \eqn{f} is different from 1 then the function is asymmetric.
}

\value{
  The value returned is a list with the following components
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
  Finney, D. J. (1979) Bioassay and the Practise of Statistical Inference, 
  \emph{Int. Statist. Rev.}, \bold{47}, 1--12.

  Seber, G. A. F. and Wild, C. J (1989) \emph{Nonlinear Regression}, New York: Wiley \& Sons (p. 330).
}

\author{Christian Ritz}

\note{
  The functions are for use with the function \code{\link{multdrc}}.
}

\seealso{
  For convenience several special cases of the function 'logistic' are available: 
  \code{\link{l2}}, \code{\link{l3}}, \code{\link{l4}} and \code{\link{l5}}. 
}

\examples{

## A four-parameter logistic model with b fixed at 1
model1 <- multdrc(FA, fct=logistic(fixed=c(1, NA, NA, NA, 1)))
summary(model1)
## Can also be specified using the function 'l4(fixed=c(1,NA,NA,NA))'

## A four-parameter logistic model with the constraint b>3
model2 <- multdrc(FA, fct=logistic(fixed=c(NA, NA, NA, NA, 1), 
lowerc=c(3, -Inf, -Inf, -Inf, -Inf)), control=mdControl(constr=TRUE))
summary(model2)

rm(model1, model2)

}
\keyword{models}
\keyword{nonlinear}
