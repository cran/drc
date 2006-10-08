\name{mlogistic}

\alias{mlogistic}
\alias{umlogistic}

%- Also NEED an '\alias' for EACH other topic documented here.

\title{The modified logistic function}

\description{
  For 'mlogistic' provides a very general way of specifying Cedergreen-Ritz-Streibig's 
  modified logistic dose response functions, under various constraints on the parameters.

  For u-shaped hormesis data: 'umlogistic' provides a very general way of specifying Cedergreen-Ritz-Streibig's 
  modified logistic dose response functions, under various constraints on the parameters.
}

\usage{
  mlogistic(lowerc = c(-Inf, -Inf, -Inf, -Inf, -Inf), 
  upperc = c(Inf, Inf, Inf, Inf, Inf), fixed = c(NA, NA, NA, NA, NA), 
  names = c("b", "c", "d", "e", "f"), alpha, scaleDose = TRUE, useDer = FALSE)
  
  umlogistic(lowerc = c(-Inf, -Inf, -Inf, -Inf, -Inf), 
  upperc = c(Inf, Inf, Inf, Inf, Inf), fixed = c(NA, NA, NA, NA, NA), 
  names = c("b", "c", "d", "e", "f"), alpha, scaleDose = TRUE, useDer = FALSE)
}
\arguments{
  \item{lowerc}{numeric vector. The lower bound on parameters. Default is minus infinity.}
  \item{upperc}{numeric vector. The upper bound on parameters. Default is plus infinity.}
  \item{fixed}{numeric vector. Specifies which parameters are fixed and at what value they are fixed. NAs for parameter that are not fixed.}
  \item{names}{a vector of character strings giving the names of the parameters (should not contain ":"). The default is reasonable (see under 'Usage').
               The order of the parameters is: b, c, d, e, f (see under 'Details').}
  \item{alpha}{numeric. The degree of hormesis. Needs to be specified!}
  \item{scaleDose}{logical. If TRUE dose values are scaled around 1 during estimation; 
                   this is required for datasets where all dose values are small. }  
  \item{useDer}{logical. If TRUE derivatives are supplied, otherwise they are not supplied. Not yet implemented!}
}
\details{
The model is given by the expression
  
  \deqn{ f(x) = c + \frac{d-c+f exp(-1/(x^{\alpha}))}{1+exp(b(log(x)-log(e)))}}
  
which is a five-parameter model (alpha is fixed).
  
It is a modification of the four-parameter logistic curve to take hormesis into account.
  
The u-shaped model is given by the expression
  
  \deqn{ f(x) = c + d - \frac{d-c+f \exp(-1/x^{\alpha})}{1+\exp(b(\log(x)-\log(e)))}}
  
}
\value{
  The value returned by the 'mlogistic' is a list with the following components
  \item{fct}{The dose response function.}
  \item{ssfct}{The self starter function.}
  \item{deriv1}{The first derivative.}  
  \item{deriv2}{The second derivative.}  
  \item{lowerc}{The lower bounds on the parameters.}
  \item{upperc}{The upper bounds on the parameters.}
  \item{edfct}{The ED function.}
  \item{sifct}{The SI function.}
  \item{maxfct}{The function for calculating the mean maximum of the dose response curve.}
}
\references{
  Cedergreen, N. and Ritz, C. and Streibig, J. C. (2005) Improved empirical models describing hormesis, 
  \emph{Environmental Toxicology and Chemistry} {\bold 24}, 3166--3172.
}

\author{Christian Ritz}

\note{This function is for use with the function \code{\link{multdrc}}.}

\seealso{Special cases of the function 'mlogistic' are \code{\link{ml3a}}, \code{\link{ml3b}}, \code{\link{ml3c}}, 
         \code{\link{ml4a}}, \code{\link{ml4b}} and \code{\link{ml4c}} where a,b and c denotes the pre-specified 
         alpha values 1, 0.5 and 0.25, respectively}

\examples{

## Modified logistic model with the constraint f>0
model1 <- multdrc(hormesis[,c(2,1)], 
fct=mlogistic(fixed=c(NA, NA, NA, NA, NA), 
lowerc=c(-Inf, -Inf, -Inf, -Inf, 0), alpha=1), 
control=mdControl(constr=TRUE))

summary(model1)
ED(model1, c(10, 50, 90))

rm(model1)


}

\keyword{models}
\keyword{nonlinear}

\concept{hormesis hormetic effect initial stimulation}
