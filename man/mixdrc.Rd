\name{mixdrc}

\alias{mixdrc}

\title{Fitting multiple dose-response curves with random effects}

\description{
  'mixdrc' fits non-linear mixed models where random effects are assigned to some or all of the parameters.
}

\usage{
mixdrc(object, random, data, startVal)
}

\arguments{
  \item{object}{an object of class 'drc'.}
  \item{random}{a character string specifying the random effects. Specification is as in \code{\link[nlme]{nlme}}.}
  \item{data}{a data frame containing the variables in the model (not optional).}
  \item{startVal}{a numeric vector of initial parameter values for the fixed effects (optional).}
}

\details{
  The function searches through a range of initial values in order to obtain convergence of the estimation procedure.
}

\value{
  An object of classes 'drc' and 'nlme'. 
}

% \references{ ~put references to the literature/web site here ~ }

\author{Christian Ritz}

%\note{It is only implemented for 'l3' and 'l4' models.}

\seealso{Non-linear mixed models can be fitted using \code{\link[nlme]{nlme}}.}

\examples{

## First fitting a model without random effects

model1 <- multdrc(SLOPE~DOSE, CURVE, data=PestSci,
collapse=data.frame(HERBICIDE, 1, 1, HERBICIDE))

## Then fitting the same model (the same fixed effects)
##  but with random effects (assigned to d)

model2 <- mixdrc(model1, random="d~1|CURVE", data=PestSci)
compParm(model2, "e", "-")

## Fitting a model without random effects
##  but with a Box-Cox transformation

model3 <- multdrc(SLOPE~DOSE, CURVE, data=PestSci,
collapse=data.frame(HERBICIDE, 1, 1, HERBICIDE), boxcox=TRUE)

## Then fitting the same model with random effects assigned to d

model4 <- mixdrc(model3, random="d~1|CURVE", data=PestSci)
compParm(model4, "e", "-")

## Fitting a model without random effects, but fixing the c parameter at 0
model5 <- multdrc(SLOPE~DOSE, CURVE, data=PestSci,
collapse=data.frame(HERBICIDE, 1, 1, HERBICIDE), boxcox=TRUE,
fct=l4(fixed=c(NA,0,NA,NA)))

## Fitting the corresponding mixed model with random effects on d
model6 <- mixdrc(model5, random="d~1|CURVE", data=PestSci)

}
\keyword{models}
\keyword{nonlinear}
