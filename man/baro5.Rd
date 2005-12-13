\name{baro5}

\alias{baro5}

\title{The modified baro5 function}

\description{
  'baro5' allows specification of the baroreflex 5-parameter dose response function,
  under various constraints on the parameters.
}

\usage{
baro5(lowerc = c(-Inf, -Inf, -Inf, -Inf, -Inf), upperc = c(Inf, Inf, Inf, Inf, Inf), fixed = c(NA, NA, NA, NA, NA), names = c("b1", "b2", "c", "d", "e"), scaleDose = TRUE, useDer = FALSE)
}

\arguments{
  \item{lowerc}{numeric vector. The lower bound on parameters. Default is minus infinity.}
  \item{upperc}{numeric vector. The upper bound on parameters. Default is plus infinity.}
  \item{fixed}{numeric vector. Specifies which parameters are fixed and at what value they are fixed. NAs for parameter that are not fixed.}
  \item{names}{a vector of character strings giving the names of the parameters (should not contain ":"). The default is reasonable (see under 'Usage').
               The order of the parameters is: b1, b2, c, d, e (see under 'Details').}
  \item{scaleDose}{logical. If TRUE dose values are scaled around 1 during estimation; 
                   this is required for datasets where all dose values are small. }
  \item{useDer}{logical. If TRUE derivatives are supplied, otherwise they are not supplied. Not implemented.}}

\details{
  The five-parameter function given by the expression

  \deqn{ y = c + \frac{d-c}{1+f\exp(b1(\log(x)-\log(e))) + (1-f)\exp(b2(\log(x)-\log(e)))}}

  \deqn{ f = 1/( 1 + \exp((2b1b2/|b1+b2|)(\log(x)-\log(e))))}
  
  If the difference between the parameters b1 and b2 is different from 0 then the function is asymmetric.
}

\value{
  The value returned by the 'baro5' is a list with the following components
  \item{fct}{The dose response function.}
  \item{ssfct}{The self starter function.}
  \item{deriv1}{The first derivative. Not available.}  
  \item{deriv2}{The second derivative. Not available.}  
  \item{lowerc}{The lower bounds on the parameters.}
  \item{upperc}{The upper bounds on the parameters.}
  \item{edfct}{The ED function. Not available.}
  \item{sifct}{The SI function. Not available.}    
}

\references{Ricketts, J. H. and Head, G. A. (1999) A five-parameter logistic equation for investigating asymmetry of
curvature in baroreflex studies. \emph{Am. J. Physiol. (Regulatory
Integrative Comp. Physiol. 46)}, \bold{277}, 441--454.}

\author{Christian Ritz}

%\note{ ~~further notes~~ }

\examples{

## See the example for the heartrate dataset

}
\keyword{models}
\keyword{nonlinear}
