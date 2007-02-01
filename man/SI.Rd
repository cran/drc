\name{SI}

\alias{SI}

\title{Comparing selectivity indices across curves}

\description{
  'SI' compares selectivity indices for arbitrary dosage across curves.
  The selectivity is the ratio between effective dosages from different curves. 
}

\usage{
  SI(object, percVec, compMatch = NULL, od = FALSE, reverse = FALSE, 
  ci = c("none", "delta", "fieller", "fls"), 
  level = ifelse(!(ci=="none"), 0.95, NULL), 
  type = c("relative", "absolute"), logBase = NULL, ...)
}

\arguments{
  \item{object}{an object of class 'drc'.}
  \item{percVec}{a numeric vector of dosage values.}
  \item{compMatch}{an optional character vector of names of assays to be compared. If not specified all comparisons are supplied.}
  \item{od}{logical. If TRUE adjustment for over-dispersion is used.}
  \item{reverse}{logical. If TRUE the order of comparison of two curves is reversed.}
  \item{ci}{character string specifying the type of confidence intervals to be supplied. The default is "none". 
  Use "delta" for asymptotics-based confidence intervals (using the delta method and the t-distribution).
  Use "fieller" for confidence intervals based on Fieller's theorem (with help from the delta method).
  Use "fls" for from logarithm scale based confidence intervals (in case the parameter in the model is 
  log(ED50) as for the \code{\link{boltzmann}}) models.}
  \item{level}{numeric. The level for the confidence intervals. Default is 0.95.}
  \item{type}{character string specifying whether absolute or relative response levels are supplied.}
  \item{logBase}{numeric. The base of the logarithm in case logarithm transformed dose values are used.}
  \item{...}{additional arguments to the function doing the calculations. 
  For instance Upper limit for bisection method. 
  Need to be larger than ED values to be calculated in ratio.}  
}

\details{  
  Fieller's theorem is incorporated using the formulas Kotz and Johnson (1983) and Finney (1978).

  For objects of class 'braincousens' or 'mlogistic' the additional argument may be the 'upper' argument
  or the 'interval' argument. The 'upper' argument specifies the upper limit of the bisection method. 
  The upper limits needs to be larger than the EDx level to be calculated. The default limit is 1000.
  The 'interval' argument should specify a rough interval in which the dose
  yielding the maximum hormetical response lies. The default interval is 'c(0.001, 1000)'. 
  Notice that the lower limit should not be set to 0 (use something like 1e-3, 1e-6, ...).
}

\value{
  A matrix with columns containing the estimates, estimated standard errors, 
  t-statistics for testing indices equal to 1 and the corresponding p-values.
}

\references{
  Finney, D. J. (1978) \emph{Statistical method in Biological Assay}, London: Charles Griffin House, 
  3rd edition (pp. 80--82).

  Kotz, S. and Johnson, N. L. (1983) \emph{Encyclopedia of Statistical Sciences Volume 3}, 
  New York: Wiley \& Sons (pp. 86--87).
}

\author{Christian Ritz}

\note{
  This function is only implemented for the following built-in functions in \code{\link{drc}}: 
  \code{\link{braincousens}}, \code{\link{cedergreen}}, \code{\link{ucedergreen}}, \code{\link{llogistic}}
  and \code{\link{weibull1}}.
}

\seealso{
  A related function is \code{\link{ED}} for calculating effect dosages.
}

\examples{

m1 <- multdrc(SLOPE~DOSE, CURVE, data=PestSci)

SI(m1, c(50,50))
SI(m1, c(10,50))
SI(m1, c(10,50), reverse = TRUE)


## Relative potency of two herbicides
m2 <- multdrc(DryMatter~Dose, Herbicide, 
data = christensenetal2003, fct = l3())

SI(m2, c(50, 50))
SI(m2, c(50, 50), ci = "fieller")

## Comparison based on absolute
##  response level

m3 <- drm(SLOPE~DOSE, CURVE,
data = PestSci, fct = LL.4())

SI(m3, c(0.5,0.5), c(2,4), type="a", ci="fieller")

SI(m3, c(55,80), c(2,4))  
# same comparison using a relative response level


## Relative potency transformed from log scale
m4 <- drm(drymat~log(dose), curve, data=cabanne[-c(1:40), ], 
pmodels = data.frame(curve,curve,1,curve), fct = B.4())

SI(m4, c(50,50), ci="fls", logBase=exp(1))

rm(m1, m2, m3, m4)

}

\keyword{models}
\keyword{nonlinear}
