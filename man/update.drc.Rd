\name{update.drc}

\alias{update.drc}

\title{Updating and re-fitting a model}

\description{
  'update' updates and re-fits a model on the basis of an object of class 'drc'.
}

\usage{

  \method{update}{drc}(object, ..., evaluate = TRUE)

}

\arguments{
  \item{object}{an object of class 'drc'.}
  \item{...}{arguments to alter in object.}
  \item{evaluate}{logical. If TRUE model is re-fit; otherwise an unevaluated call is returned.}
}

\value{
  An object of class 'drc'.
}

\author{Christian Ritz}

\examples{

## Fitting a four-parameter Weibull model
model1 <- multdrc(FA, fct = w4())

## Updating 'model1' by fitting a three-parameter Weibull model
model2 <- update(model1, fct = w3())
anova(model2, model1)


rm(model1, model2)
}
\keyword{models}
\keyword{nonlinear}
