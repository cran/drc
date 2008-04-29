\name{comped}

\alias{comped}

\title{Comparison of effective dose values}

\description{
  Comparison of a pair of effective dose values from independent experiments 
  where only the estimates and their standard errors are reported.
}

\usage{
  comped(est, se, log = TRUE, interval = TRUE, operator = c("-", "/"), level = 0.95, df = NULL)
}

\arguments{
  \item{est}{a numeric vector of length 2 containing the two estimates}
  \item{se}{a numeric vector of length 2 containing the two standard errors}
  \item{log}{logical indicating whether or not estimates and standard errors are on log scale}
  \item{interval}{logical indicating whether or not a confidence interval should be returned}
  \item{operator}{character string taking one of the two values "-" (default) or "/" corresponding to a comparison
  based on the difference or the ratio.}
  \item{level}{numeric value giving the confidence level}
  \item{df}{numeric value specifying the degrees of freedom for the percentile used in the confidence interval (optional)}
}

\details{
   The choice "/" for the argument \code{operator} and FALSE for \code{log} will result in estimation of a socalled
   relative potency (sometimes also called a selectivity index).

   The combination TRUE for \code{log} and "/" for \code{operator} only influences the confidence interval, 
   that is no ratio is calculated based on logarithm-transformed effective dose values.
   
   By default confidence interval relies on percentiles in the normal distribution.  
   
   In case the entire dataset is available the functions \code{\link{drm}} and (subsequently) \code{\link{SI}}
   should be used.
}

\value{
  A matrix with the estimated difference or ratio and the associated standard error and the resulting confidence
  interval (unless not requested).
}

\references{
  Wheeler, M. W. and Park, R. M. and Bailer, A. J. (2006)
  Comparing median lethal concentration values using confidence interval overlap or ratio tests,
  \emph{Environmental Toxicology and Chemistry}, \bold{25}, 1441--1441.
}

\author{Christian Ritz}

\note{
   The development of the function \code{comped} is offspin of a project on statistical analysis of toxicity data 
   funded by the Danish EPA ("Statistisk analyse og biologisk tolkning af toksicitetsdata", MST j.nr. 669-00079).
}

\seealso{
  The function \code{\link{ED}} calculates arbitrary effective dose values based on a model fit. The function
  \code{\link{SI}} calculates relative potencies based on arbitrary effective dose values.
}

\examples{

comped(c(21.173, 44.718), sqrt(c(11.87, 8.42)), log=FALSE, operator="/", level=0.99)

comped(c(21.173, 44.718), sqrt(c(11.87, 8.42)), log=FALSE, operator="/", level=0.99, interval=FALSE)

comped(c(21.173, 44.718), sqrt(c(11.87, 8.42)), log=FALSE, operator="-", level=0.99)

}
\keyword{models}
\keyword{nonlinear}
