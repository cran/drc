\name{SI}

\alias{SI}
\alias{relpot}

\title{Comparing selectivity indices across curves}

\description{
  'SI' compares selectivity indices for arbitrary dosage across curves.
  The selectivity is the ratio between effective dosages from different curves. 
}

\usage{
  SI(object, percVec, compMatch = NULL, od = FALSE, reverse = FALSE, 
  ci = c("none", "delta", "fieller", "fls"), 
  level = ifelse(!(ci=="none"), 0.95, NULL),
  reference = c("control", "upper"),  
  type = c("relative", "absolute"), logBase = NULL,
  display = TRUE, ...)
  
  relpot(object, plotit = TRUE, compMatch = NULL, percVec = NULL, ci = "none", 
  type = c("relative", "absolute"), scale = c("original", "percent", "unconstrained"), ...)
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
  Use "fls" for from-logarithm-scale-based confidence intervals (in case the parameter in the model fit is 
  log(ED50) as is the case for the \code{\link{boltzmann}} or \code{\link{llogistic2}} models); currently
  the argumetn \code{logBase} then also needs to be specified.}
  \item{level}{numeric. The level for the confidence intervals. Default is 0.95.}
  \item{reference}{character string. Is the upper limit or the control level the reference?}  
  \item{type}{character string specifying whether absolute or relative response levels are supplied.}
  \item{logBase}{numeric. The base of the logarithm in case logarithm transformed dose values are used.}
  \item{display}{logical. If TRUE results are displayed. Otherwise they are not (useful in simulations).}  
  \item{...}{In \code{SI}: additional arguments to the function doing the calculations. 
  For instance the upper limit for the bisection method 
  needs to be larger than the ED values used in the required relative pontency.
  In \code{relpot}: additional graphical parameters.}  
  \item{plotit}{logical. If TRUE the relative potencies are plotted as a function of the response level.}    
  \item{scale}{character string indicating the scale to be used on the x axis: original or percent response level
  (only having an effect for type="relative").}    
}

\details{
  The function \code{relpot} is a convenience function, which is useful for assessing how the relative potency
  changes as a function of the response level (eg for plotting as outlined in Ritz \emph{et al} (2006)).
  
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
  
  Ritz, C. and Cedergreen, N. and Jensen, J. E. and Streibig, J. C. (2006) 
  Relative potency in nonsimilar dose-response curves, \emph{Weed Science}, \bold{54}, 407--412.
}

\author{Christian Ritz}

\note{
  This function is only implemented for the following built-in functions available in the package \emph{drc}: 
  \code{\link{braincousens}}, \code{\link{cedergreen}}, \code{\link{ucedergreen}}, \code{\link{llogistic}}
  and \code{\link{weibull1}}.
}

\seealso{
  A related functions \code{\link{ED}} (for calculating effective doses) and 
  \code{\link{relpot}} (for displaying relative potencies).
}

\examples{

m1 <- drm(SLOPE~DOSE, CURVE, data = spinach, fct = LL.4())

SI(m1, c(50,50))
SI(m1, c(10,50))
SI(m1, c(10,50), reverse = TRUE)


## Relative potency of two herbicides
m2 <- drm(DryMatter~Dose, Herbicide, 
data = S.alba, fct = LL.3())

SI(m2, c(50, 50))  
SI(m2, c(50, 50), ci = "delta")
SI(m2, c(50, 50), ci = "fieller")

## Comparison based on absolute
##  response level

m3 <- drm(SLOPE~DOSE, CURVE,
data = spinach, fct = LL.4())

SI(m3, c(0.5,0.5), c(2,4), type="a", ci="fieller")

SI(m3, c(55,80), c(2,4))  
# same comparison using a relative response level


## Relative potency transformed from log scale
m4 <- drm(drymatter~log(dose), treatment, data=G.aparine[-c(1:40), ], 
pmodels = data.frame(treatment,treatment,1,treatment), fct = LL2.4())

SI(m4, c(50,50), ci="fls", logBase=exp(1))

}

\keyword{models}
\keyword{nonlinear}
