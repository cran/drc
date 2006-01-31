\name{mixdrc}
\alias{mixdrc}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{Fitting multiple dose-response curves with random effects}

\description{
  'mixdrc' fits non-linear mixed models where random effects are assigned to some or all of the parameters.
}

\usage{
mixdrc(object, random, data)
}

\arguments{
  \item{object}{an object of class 'drc'.}
  \item{random}{a character string specifying the random effects. Specification is as in \code{\link{nlme}}.}
  \item{data}{a data frame containing the variables in the model (not optional).}
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

%\seealso{ ~~objects to See Also as \code{\link{~~fun~~}}, ~~~ }

\examples{

## First fitting model without random effects
model1 <- multdrc(SLOPE~DOSE, CURVE, data=PestSci,
collapse=data.frame(HERBICIDE, 1, 1, HERBICIDE))

## Then fitting model with random effects (assigned to d)
model2 <- mixdrc(model1, random="d~1|CURVE", data=PestSci)
compParm(model2, "e", "-")

## First fitting model without random effects
##  but with a Box-Cox transformation
model3 <- multdrc(SLOPE~DOSE, CURVE, data=PestSci,
collapse=data.frame(HERBICIDE, 1, 1, HERBICIDE), boxcox=TRUE)

## Then fitting model with random effects (assigned to d)
model4 <- mixdrc(model3, random="d~1|CURVE", data=PestSci)
compParm(model4, "e", "-")


rm(model1, model2, model3, model4)

}
\keyword{models}
\keyword{nonlinear}
