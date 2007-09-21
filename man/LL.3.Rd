\name{LL.3}

\alias{LL.3}
\alias{LL.3u}

\alias{l3}
\alias{l3u}

\alias{LL2.3}
\alias{LL2.3u}

\title{The three-parameter log-logistic function}

\description{
  'LL.3' and 'LL2.3' provide the three-parameter log-logistic function where the lower limit is equal to 0. 
  
  'LL.3u' and 'LL2.3u' provide three-parameter logistic function where the upper limit is equal to 1, mainly
  for use with binomial/quantal response. 
}

\usage{  
  LL.3(fixed = c(NA, NA, NA), names = c("b", "d", "e"), ...)
  
  LL.3u(fixed = c(NA, NA, NA), names = c("b", "c", "e"), ...)
  
  l3(fixed = c(NA, NA, NA), names = c("b", "d", "e"), ...)

  l3u(fixed = c(NA, NA, NA), names = c("b", "c", "e"), ...)  
  
  LL2.3(fixed = c(NA, NA, NA), names = c("b", "d", "e"), ...)
  
  LL2.3u(fixed = c(NA, NA, NA), names = c("b", "c", "e"), ...)  
}

\arguments{
  \item{fixed}{numeric vector. Specifies which parameters are fixed and at what value they are fixed. 
  NAs for parameter that are not fixed.}
  \item{names}{a vector of character strings giving the names of the parameters. The default is reasonable.}
  \item{...}{Additional arguments (see \code{\link{llogistic}}).}
}

\details{
  The three-parameter logistic function with lower limit 0 is
  \deqn{ f(x) = 0 + \frac{d-0}{1+\exp(b(\log(x)-\log(e)))}}
  
  or in another parameterisation
  \deqn{ f(x) = 0 + \frac{d-0}{1+\exp(b(\log(x)-e))}}
    
  The three-parameter logistic function with upper limit 1 is
  \deqn{ f(x) = c + \frac{1-c}{1+\exp(b(\log(x)-\log(e)))}}
  
  or in another parameterisation
  \deqn{ f(x) = c + \frac{1-c}{1+\exp(b(\log(x)-e))}}
  
  Both functions are symmetric about the inflection point (\eqn{e}).  
}

\value{
  See \code{\link{llogistic}}.
}

\references{
  Finney, D. J. (1971) \emph{Probit Analysis}, Cambridge: Cambridge University Press.
}

\author{Christian Ritz}

\note{
  This function is for use with the function \code{\link{multdrc}}.
}

\seealso{
  Related functions are \code{\link{LL.2}}, \code{\link{LL.4}}, \code{\link{LL.5}} and the more general 
  \code{\link{llogistic}}.
}

\examples{

## Fitting model with lower limit equal 0
model1 <- multdrc(ryegrass, fct=LL.3())
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
logDose=10, fct=LL.3u(), type="binomial", 
collapse=data.frame(treatment, 1, treatment))

summary(model2)
anova(model2)
plot(model2, conLevel=-1, ylim=c(0.1, 1.3))
abline(h=1, lty=2)

rm(model1, model2, n, r, treatment, finney_ex19)

}
\keyword{models}
\keyword{nonlinear}
