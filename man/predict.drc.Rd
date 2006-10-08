\name{predict.drc}

\alias{predict.drc}

\title{Prediction}

\description{
  Predicted values for models of class 'drc'.
}

\usage{

  \method{predict}{drc}(object, newdata, se.fit = TRUE, interval = c("none", "confidence", "prediction"), 
  level = 0.95, na.action = na.pass, od = FALSE, ...)

}

\arguments{
  \item{object}{an object of class 'drc'.}
  \item{newdata}{An optional data frame in which to look for variables with which to predict.
  If omitted, the fitted values are used.}
  \item{se.fit}{logical. If TRUE standard errors are required.}
  \item{interval}{character string. Type of interval calculation: "none", "confidence" or "prediction".}
  \item{level}{Tolerance/confidence level.}
  \item{na.action}{function determining what should be done with missing values in 'newdata'. 
  The default is to predict 'NA'.}
  \item{od}{logical. If TRUE adjustment for over-dispersion is used.}
  \item{\dots}{further arguments passed to or from other methods.}
}

\details{
  For the built-in log-logistics and Weibull-type models standard errors and confidence/prediction intervals 
  can be calculated. At the moment it only works for the situations where all observations are assumed to 
  have a common variance.
}

\value{
  A matrix with as many rows as there are dose values provided in 'newdata' or in the original dataset 
  (in case 'newdata' is not specified) and columns with fitted, standard errors, lower and upper limits
  of confidence intervals.
}

%\references{}

\author{Christian Ritz}

%\note{}

\seealso{For details are found in the help page for \code{\link{predict.lm}}.}

\examples{

## Fitting a model
model1 <- multdrc(SLOPE~DOSE, CURVE, data=PestSci)

## Predicting values a dose=2 (with standard errors)
predict(model1, data.frame(dose=2, CURVE=c("1", "2", "3")))

## Getting confidence intervals
predict(model1, data.frame(dose=2, CURVE=c("1", "2", "3")), 
interval = "confidence")

## Getting prediction intervals
predict(model1, data.frame(dose=2, CURVE=c("1", "2", "3")), 
interval = "prediction")

rm(model1)

}

\keyword{models}
\keyword{nonlinear}
