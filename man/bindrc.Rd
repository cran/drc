\name{bindrc}
\alias{bindrc}


\title{Fitting multiple dose-response curves for binomial data}

\description{
  'bindrc' fits logit/probit regression models to multiple dose-response curves.
  The models are linear in the logarithm of dose in the logit/probit scale.
}

\usage{
bindrc(formula, weights, curve, int, slope, data = NULL, na.action = na.fail, 
link = "logit", log =TRUE, startVal, lower, upper, fixedLow = NULL, fixedUp = NULL, 
names = c("a", "b", "c", "d"))
}


\arguments{
  \item{formula}{a symbolic description of the model to be fit in the form 'response \eqn{~} dose'.}
  \item{weights}{a numeric containing number of trials per experiment, eg per dose level.}
  \item{curve}{a numeric vector or factor containing the grouping of the data into separate curves.
               If not specified all observations are assumed to belong to a single curve!}
  \item{int}{a numeric vector or factor containing the grouping of the data for the intercept.}
  \item{slope}{a numeric vector or factor containing the grouping of the data for the intercept.}
  \item{data}{an optional data frame containing the variables in the model.}
  \item{na.action}{a function which indicates what should happen when the data contain 'NA's.  The default is 'na.fail'. To omit 'NA's use
  'na.omit'.}
  \item{link}{a character string specifying the model: "logit", "probit" or "cloglog".}
  \item{log}{logical. If TRUE (default) independent variable is log-transformed. If FALSE no log-transformation is used.}
  \item{startVal}{an optional numeric vector containing start values for all parameters in the model. 
                  Overrules any self starter facility.}
  \item{lower}{a numeric vector or factor containing the grouping of the data for the lower limit.}
  \item{upper}{a numeric vector or factor containing the grouping of the data for the upper limit.}
  \item{fixedLow}{a vector for fixing the lower limit for some curves.}
  \item{fixedUp}{a vector for fixing the upper limit for some curves.}
  \item{names}{a character vector with the names of the parameters.}
}

\details{
  The model fit has the form
  \deqn{ p = c + (d-c) F( \alpha + \beta \log dose)}
  for each curve with F the inverse of logit or probit.
  
  Maximum likelihood estimation is employed.
}

\value{
  An object of class 'drc'. Various extractors exist for objects of this class.
}

%\references{}

\author{Christian Ritz}

%\note{}

\seealso{See \code{\link{multdrc}} for fitting data with continuous response.}

\examples{

}
\keyword{models}
\keyword{nonlinear}
