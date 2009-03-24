\name{BC.5}

\alias{BC.5}
\alias{bcl4}

\alias{BC.4}
\alias{bcl3}

\title{The Brain-Cousens hormesis models}

\description{
  'BC.4' and 'BC.5' provide the Brain-Cousens modified log-logistic models for describing u-shaped hormesis.
}

\usage{
  BC.5(fixed = c(NA, NA, NA, NA, NA), names = c("b", "c", "d", "e", "f"))

  BC.4(fixed = c(NA, NA, NA, NA), names = c("b", "d", "e", "f"))
}

\arguments{
  \item{fixed}{numeric vector specifying which parameters are fixed and at which values they are fixed. 
  NAs designate parameters that are not fixed.}
  \item{names}{a vector of character strings giving the names of the parameters. The default is reasonable.}
}

\details{
  The Brain-Cousens model is given by the expression
  \deqn{ f(x) = c + \frac{d-c+fx}{1+\exp(b(\log(x)-\log(e)))}}
  which is a five-parameter model.
  
  It is a modification of the log-logistic curve to take u-shaped hormesis into account.
  
  Fixing the lower limit at 0 yields the four-parameter model 
  \deqn{ f(x) = 0 + \frac{d-0+fx}{1+\exp(b(\log(x)-\log(e)))}}
}

\value{
  See \code{\link{braincousens}}.
}

\references{
  Brain, P. and Cousens, R. (1989) An equation to describe dose responses 
  where there is stimulation of growth at low doses, 
  \emph{Weed Research}, \bold{29}, 93--96.

  van Ewijk, P. H. and Hoekstra, J. A. (1993) 
  Calculation of the EC50 and its Confidence Interval When Subtoxic Stimulus Is Present,
  \emph{Ecotoxicology and Environmental Safety}, \bold{25}, 25--32.  
}

\author{Christian Ritz}

\note{This function is for use with the function \code{\link{drm}}.}

\seealso{
  More details are found for the general model function \code{\link{braincousens}}.
}

\examples{

## Fitting the data in van Ewijk and Hoekstra (1993)
lettuce.bcm1 <- drm(weight ~ conc, data = lettuce, fct=BC.5())
modelFit(lettuce.bcm1)
plot(lettuce.bcm1)

lettuce.bcm2 <- drm(weight ~conc, data = lettuce, fct=BC.4())
summary(lettuce.bcm2)
ED(lettuce.bcm2, c(50))  
# compare the parameter estimate and 
# its estimated standard error 
# to the values in the paper by 
# van Ewijk and Hoekstra (1993)


## Brain-Cousens model with the constraint b>3
ryegrass.bcm1 <- drm(rootl ~conc, data = ryegrass, fct = BC.5(), 
lower = c(3, -Inf, -Inf, -Inf, -Inf), control = drmc(constr=TRUE))

summary(ryegrass.bcm1)

## Brain-Cousens model with the constraint f>0 
## (no effect as the estimate of f is positive anyway)
ryegrass.bcm2 <- drm(rootl ~conc, data = ryegrass, fct = BC.5(), 
lower = c(-Inf, -Inf, -Inf, -Inf, 0), control = drmc(constr=TRUE))

summary(ryegrass.bcm2)

}
\keyword{models}
\keyword{nonlinear}
