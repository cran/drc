\name{anova.drc}

\alias{anova.drc}

\alias{anova.drclist}

\title{ANOVA for non-linear model fits}

\description{
  'anova' produces an analysis of variance table for one or two non-linear model fits.
}

\usage{
anova.drc(object, ..., test = NULL)

anova.drclist(object, ..., test = NULL)
}

\arguments{
  \item{object}{an object of class 'drc'.}
  \item{...}{additional arguments.}
  \item{test}{a character string specifying the test statistic to be applied.
              Use "od" to assess overdispersion for binomial data.}
}

\details{
  Specifying only a single object gives a test for lack-of-fit, comparing the non-linear regression
  model to a more general one-way or two-way ANOVA model.
  
  If two objects are specified a test for reduction from the larger to the smaller model is given. (This only makes statistical 
  sense if the models are nested, that is: one model is a submodel of the other model.)
}

\value{
  An object of class 'anova'.
}

\references{Bates, D. M. and Watts, D. G. (1988) \emph{Nonlinear Regression Analysis and Its Applications}, New York: Wiley \& Sons (pp. 103--104).} 

\author{Christian Ritz}

\seealso{The function \code{\link{anova.lm}} for linear models.}

\examples{

## Comparing the nonlinear model to a one-way ANOVA model using an F test
model1 <- multdrc(FA, fct = w4())
anova(model1)   

## Comparing a Gompertz three- and four-parameter models using an F test
model2 <- multdrc(FA, fct = w3())
anova(model2, model1)  # reduction to 'w3' not possible (highly significant)   


rm(model1, model2)
}
\keyword{models}
\keyword{nonlinear}
