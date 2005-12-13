\name{coef.drc}
\alias{coef.drc}
\title{Extract Model Coefficients}
\description{
   Extract parameter estimates.
}
\usage{
coef.drc(object, ...)
}
\arguments{
  \item{object}{an object of class 'drc'.}
  \item{...}{additional arguments.}
}
\value{
  A vector of parameter coefficients which are extracted from the model object 'object'.
}
\author{Christian Ritz}
\note{This function may work even in cases where 'summary' does not, because the parameter estimates 
      are retrieved directly from the model fit object without any additional computations of summary statistics and standard errors.}

\seealso{A more comprehensive summary is obtained using \code{\link{summary.drc}}.}
\examples{

model1 <- multdrc(FA)
coef(model1)

model2 <- multdrc(FA, boxcox=TRUE)
summary(model2)  # not a satisfactory fit

## Usinf estimates from 'model1' as start values
model3 <- multdrc(FA, boxcox=TRUE, startVal=coef(model1))
summary(model3)  # a satisfactory fit

rm(model1, model2, model3)

}
\keyword{models}
\keyword{nonlinear}
