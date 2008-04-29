\name{MM}

\alias{MM.2}
\alias{MM.3}

\title{Michaelis-Menten model}

\description{
  The functions can be used to fit (shifted) Michaelis-Menten models that are used
  for modeling enzyme kinetics, weed densities etc.
}

\usage{ 
  MM.2(fixed = c(NA, NA), names = c("Vm", "K"))
  
  MM.3(fixed = c(NA, NA, NA), names = c("y0", "Vm", "K"))  
}

\arguments{
  \item{fixed}{numeric vector. Specifies which parameters are fixed and at what value they are fixed. 
  NAs for parameter that are not fixed.}
  \item{names}{a vector of character strings giving the names of the parameters (should not contain ":"). 
  The order of the parameters is: y0, Vm, K (see under 'Details').}}

\details{
  The model is given by the expression
  
  \deqn{ f(x) = y0 + \frac{Vm x}{K + x}}
  
  It is a three-parameter model. A commonly used derived two-parameter model (\code{MM.2}) is obtained by
  setting y0=0.
}

\value{
  A list of class \code{drcMean}, containing the mean function, the self starter function,
  the parameter names and other components such as derivatives and a function for calculating ED values.
}

%\references{}

\author{Christian Ritz}

\note{
  At the moment the implementation cannot deal with infinite concentrations.
}

\seealso{Related models are the asymptotic regression models \code{\link{AR.2}} and \code{\link{AR.3}}.}

\examples{

met.mm.m1<-drm(gain~dose, product, data=methionine, fct=MM.3(), pmodels = list(~1, ~factor(product), ~factor(product)))
plot(met.mm.m1, log = "", ylim=c(1450, 1800))
summary(met.mm.m1)

## Calculating bioefficacy: approach 1
coef(met.mm.m1)[4] / coef(met.mm.m1)[5] * 100

## Calculating bioefficacy: approach 2
SI(met.mm.m1, c(50,50))

## Simplified models
met.mm.m2<-drm(gain~dose, product, data=methionine, fct=MM.3(), pmodels = list(~1, ~factor(product), ~1))
anova(met.mm.m2, met.mm.m1)  # model reduction not possible

met.mm.m3<-drm(gain~dose, product, data=methionine, fct=MM.3(), pmodels = list(~1, ~1, ~factor(product)))
anova(met.mm.m3, met.mm.m1)  # model reduction not possible

}
\keyword{models}
\keyword{nonlinear}
