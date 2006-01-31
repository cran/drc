\name{compParm}

\alias{compParm}

\title{Comparison of parameters}

\description{
  Compare parameters from different assays, either by means of ratios or differences.
}

\usage{
compParm(object, strVal, operator = "/", od = FALSE)
}

\arguments{
  \item{object}{an object of class 'drc'.}
  \item{strVal}{a name of parameter to compare.}
  \item{operator}{a character. If equal to "/" (default) parameter ratios are compared. If equal to "-" parameter differences are compared.}
  \item{od}{logical. If TRUE adjustment for over-dispersion is used.}
}

\value{
  A matrix with columns containing the estimates, estimated standard errors, values of t-statistics and p-values for the null hypothesis that 
  the ratio equals 1 or that the difference equals 0 (depending on the \code{operator} argument). 
  
  No adjustment for multiplicity is applied. See \code{\link{p.adjust}} on how to do adjust p-values.
}

\details{
  The function compares actual parameter estimates, and therefore the results depend on the parameterisation used. Probably it is most useful
  in combination with 'collapse' arguments (in \code{\link{multdrc}}) that are data frames or lists with formulas without intercept (-1). 
}

\author{Christian Ritz}

\examples{

# Fitting a model with names assigned to the parameters!
model1 <- multdrc(SLOPE~DOSE, CURVE, data=PestSci,
fct=l4(names = c("b", "lower", "upper", "ed50")))

## Calculating ratios of parameter estimates for the parameter named "ed50"
compParm(model1, "ed50")      

## Calculating differences between parameter estimates for the parameter named "ed50"
compParm(model1, "ed50", "-")  

rm(model1)

}
\keyword{models}
\keyword{nonlinear}
