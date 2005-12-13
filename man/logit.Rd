\name{logit}
\alias{logit}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{The logit function}

\description{
  'logit' provides a very general way of specifying logit dose-response functions for binomial data.
}

\usage{
logit(lowerc = c(-Inf, -Inf, -Inf, -Inf), upperc = c(Inf, Inf, Inf, Inf), 
fixed = c(NA, NA, 0, 1), names = c("a", "b", "c", "d"), 
useDer = FALSE, logTrans = TRUE)
}

\arguments{
  \item{lowerc}{numeric vector. The lower bound on parameters. Default is minus infinity.}
  \item{upperc}{numeric vector. The upper bound on parameters. Default is plus infinity.}
  \item{fixed}{numeric vector. Specifies which parameters are fixed and at what value they are fixed. NAs for parameter that are not fixed.}
  \item{names}{a vector of character strings giving the names of the parameters (should not contain ":"). The default is reasonable (see under 'Usage').
               The order of the parameters is: b, c, d, e, f (see under 'Details').}
  \item{useDer}{logical. If TRUE derivatives are supplied, otherwise they are not supplied. Not yet implemented!}
  \item{logTrans}{logical. If TRUE (default) the logarithm of the dose is used in the linear predictor.}
}

\details{
  The default arguments yields the five-parameter logistic function given by the expression
  \deqn{ f(x) = c + (d - c) inverse logit (a + b \log(x))}

}

\value{
  The value returned by the 'logistic' is a list with the following components
  \item{fct}{The dose response function.}
  \item{ssfct}{The self starter function.}
  \item{deriv1}{The first derivative.}  
  \item{deriv2}{The second derivative.}  
  \item{lowerc}{The lower bounds on the parameters.}
  \item{upperc}{The upper bounds on the parameters.}
  \item{edfct}{The ED function.}
  \item{sifct}{The SI function (not yet available).}    
}

\references{
   Finney, D. J. (1971) \emph{Probit Analysis}, Cambridge: Cambridge University Press.
}

\author{Christian Ritz}

%\note{ ~~further notes~~ }

\seealso{Special cases of the function 'logistic' are \code{\link{l3}}, \code{\link{l4}} and \code{\link{l5}}.}

\examples{

## Example dataset from Finney (1971) - example 19
logdose <- c(2.17, 2,1.68,1.08,-Inf,1.79,1.66,1.49,1.17,0.57)
n <- c(142,127,128,126,129,125,117,127,51,132)
r <- c(142,126,115,58,21,125,115,114,40,37)
treatment <- factor(c("w213","w213","w213","w213","control","w214","w214","w214","w214","w214"))
finney_ex19 <-data.frame(logdose, n, r, treatment)

## Fitting logistic regression model
model1 <- multdrc(r/n~logdose, treatment, logDose=10, weights=n, fct=logit(), data=finney_ex19[-5,])
anova(model1)
summary(model1)

## Fitting extended logistic regression model (lower limit above 0 is estimated)
model2 <- multdrc(r/n~logdose, treatment, logDose=10, weights=n, 
fct=logit(fixed=c(NA,NA,NA,1)), cm="control", 
collapse=data.frame(treatment,treatment,1), data=finney_ex19)
anova(model2)
summary(model2)

rm(logdose, n, r, treatment, finney_ex19, model1, model2)

}
\keyword{models}
\keyword{nonlinear}
