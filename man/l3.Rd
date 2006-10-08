\name{l3}

\alias{l3}
\alias{l3u}

\title{The three-parameter logistic function}

\description{
  'l3' provides the three-parameter logistic function where the lower limit is equal to 0. 
  
  'l3u' provides three-parameter logistic function where the upper limit is equal to 1, mainly
  for use with binomial/quantal response. 
}

\usage{
  l3(fixed = c(NA, NA, NA), names = c("b", "d", "e"), useDer = FALSE, w =FALSE)

  l3u(fixed = c(NA, NA, NA), names = c("b", "c", "e"), useDer = FALSE, w = FALSE)
}

\arguments{
  \item{fixed}{numeric vector. Specifies which parameters are fixed and at what value they are fixed. 
  NAs for parameter that are not fixed.}
  \item{names}{a vector of character strings giving the names of the parameters. The default is reasonable.}
  \item{useDer}{logical. If TRUE derivatives are supplied, otherwise they are not supplied. Not implemented!}
  \item{w}{logical. If TRUE weighted regression is applied to determine the start values.}
}

\details{
  The three-parameter logistic function with lower limit 0 is
  \deqn{ f(x) = 0 + \frac{d-0}{1+\exp(b(\log(x)-\log(e)))}}
  
  The three-parameter logistic function with upper limit 1 is
  \deqn{ f(x) = c + \frac{1-c}{1+\exp(b(\log(x)-\log(e)))}}
  
  Both functions are symmetric about the inflection point (\eqn{e}).  
}

\value{
  See \code{\link{logistic}}.
}

\references{
  Finney, D. J. (1971) \emph{Probit Analysis}, Cambridge: Cambridge University Press.
}

\author{Christian Ritz}

\note{
  This function is for use with the function \code{\link{multdrc}}.
}

\seealso{
  Related functions are \code{\link{l2}}, \code{\link{l4}}, \code{\link{l5}} and the more general 
  \code{\link{logistic}}.
}

\examples{

## Fitting model with lower limit equal 0
model1 <- multdrc(FA, fct=l3())
summary(model1)

## Fitting binomial response
##  with non-zero control response

## Example dataset from Finney (1971) - example 19
logdose <- c(2.17, 2,1.68,1.08,-Inf,1.79,1.66,1.49,1.17,0.57)
n <- c(142,127,128,126,129,125,117,127,51,132)
r <- c(142,126,115,58,21,125,115,114,40,37)
treatment <- factor(c("w213","w213","w213","w213",
"control","w214","w214","w214","w214","w214"))
finney_ex19 <- data.frame(logdose, n, r, treatment)

## Fitting model where the lower limit is estimated
model2 <- multdrc(r/n~logdose, treatment, weights=n, data=finney_ex19, 
logDose=10, fct=l3u(), type="binomial", 
collapse=data.frame(treatment, 1, treatment))

summary(model2)
anova(model2)
plot(model2, conLevel=-1, ylim=c(0.1, 1.3))
abline(h=1, lty=2)

rm(model1, model2, n, r, treatment, finney_ex19)

}
\keyword{models}
\keyword{nonlinear}
