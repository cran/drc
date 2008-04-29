\name{asymreg}

\alias{AR.2}
\alias{AR.3}

\title{Asymptotic regression model}

\description{
  Providing the mean function and the corresponding self starter function for the asymptotic regression model.
}

\usage{

AR.2(fixed = c(NA, NA), 
names = c("b", "c"))

AR.3(fixed = c(NA, NA, NA), 
names = c("a", "b", "c"))

}

\arguments{
  \item{fixed}{numeric vector. Specifies which parameters are fixed and at what value they are fixed. 
  NAs for parameter that are not fixed.}
  \item{names}{vector of character strings giving the names of the parameters (should not contain ":").}
}

\details{
  The asymptotic regression model is a three-parameter model with mean function:
  
  \deqn{ f(x) = a + (b-a)(1-\exp(-c x))}
  
  The parameter a is the lower limit (at x=0), the parameter b is the upper limit 
  and the parameter c (>0) is determining the steepness of the increase as x.
}

\value{
  A list of class \code{drcMean}, containing the mean function, the self starter function,
  the parameter names and other components such as derivatives and a function for calculating ED values.
}

%\references{ ~put references to the literature/web site here ~ }

\author{Christian Ritz}

\note{
  The functions are for use with the function \code{\link{drm}}.
}

\seealso{A very similar, but decaying model is the exponential decay model \code{\link{expDecay}}.}

\examples{

## First model
met.as.m1<-drm(gain~dose, product, data=methionine, fct=AR.3(), 
pmodels = list(~1, ~factor(product), ~factor(product)))
plot(met.as.m1, log = "", ylim=c(1450, 1800))
summary(met.as.m1)

## Calculating bioefficacy: approach 1
coef(met.as.m1)[5] / coef(met.as.m1)[4] * 100

## Calculating bioefficacy: approach 2
SI(met.as.m1, c(50,50))

## Simplified models
met.as.m2<-drm(gain~dose, product, data=methionine, fct=AR.3(), 
pmodels = list(~1, ~1, ~factor(product)))
anova(met.as.m2, met.as.m1)  # simplification not possible

met.as.m3<-drm(gain~dose, product, data=methionine, fct=AR.3(), pmodels = list(~1, ~factor(product), ~1))
anova(met.as.m3, met.as.m1)  # simplification not possible

}
\keyword{models}
\keyword{nonlinear}

