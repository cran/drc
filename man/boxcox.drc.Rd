\name{boxcox.drc}

\alias{boxcox.drc}

\title{Transform-both-sides Box-Cox transformation}

\description{
  Finds the optimal Box-Cox transformation for non-linear regression.
}

\usage{
boxcox.drc(object, lambda = seq(-2, 2, 1/10), plotit = TRUE, eps = 1/50, bcAdd = 0, level = 0.95, xlab = expression(lambda), ylab = "log-Likelihood", ...)
}

\arguments{
  \item{object}{object of class \code{drc}.}
  \item{lambda}{numeric vector of lambda values; the default is (-2, 2) in steps of 0.1.}
  \item{plotit}{logical which controls whether the result should be plotted.}
  \item{eps}{numeric value: the tolerance for lambda = 0; defaults to 0.02.}
  \item{bcAdd}{numeric value specifying the constant to be added on both sides prior to Box-Cox transformation. 
  The default is 0.}
  \item{level}{numeric value: the confidence level required.}
  \item{xlab}{character string: the label on the x axis, defaults to "lambda".}
  \item{ylab}{character string: the label on the y axis, defaults to "log-likelihood".}
  \item{\dots}{additional graphical parameters.}
}

\details{
  The optimal lambda value is determined using a profile likelihood approach: 
  For each lambda value the non-linear regression model is fitted and the lambda value resulting in thre largest
  value of the log likelihood function is picked.
}

\value{
  An object of class \code{\link{nls}} (returned invisibly). 
  If plotit = TRUE a plot of loglik vs lambda is shown indicating a confidence interval (by default 95%) about 
  the optimal lambda value.
}

\references{
    Carroll, R. J. and Ruppert, D. (1988) \emph{Transformation and Weighting in Regression}, 
    New York: Chapman and Hall (Chapter 4).
}

\author{Christian Ritz}

\note{}

\seealso{
  For linear regression the analogue is \code{\link[MASS]{boxcox}}. 
  
%  For general non-linear regression models
%  fitted using \code{\link{nls}} you can use \code{\link{boxcox.nls}}.
}

\examples{

## Fitting log-logistic model without transformation
m1 <- drm(ryegrass, fct = LL.4())
summary(m1)

## Fitting the same model with optimal Box-Cox transformation
m2 <- boxcox(m1)
summary(m2)

}

\keyword{models}
\keyword{nonlinear}
