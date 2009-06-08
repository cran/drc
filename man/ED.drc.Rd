\name{ED.drc}

\alias{ED}
\alias{ED.drc}
\alias{ED.mrdrc}

\title{Estimating effective doses}

\description{

  \code{ED} estimates effective doses (ECp/EDp/ICp) for given reponse levels.
}

\usage{
  \method{ED}{drc}(object, respLev, interval = c("none", "delta", "fls", "tfls"), 
  level = ifelse(!(interval == "none"), 0.95, NULL),
  reference = c("control", "upper"), type = c("relative", "absolute"), lref, uref,
  bound = TRUE, od = FALSE, display = TRUE, pool = TRUE, logBase = NULL, ...)
  
  \method{ED}{mrdrc}(object, respLev, interval = c("none", "approximate", "bootstrap"), level = 0.95, 
  method = c("bisection", "grid"), cgridsize = 20, gridsize = 100, display = TRUE, lower, upper,
  intType = c("confidence", "prediction"), minmax = c("response", "dose"), n = 1000, seedVal = 200810311, ...) 
}

\arguments{
  \item{object}{an object of class 'drc'.}
  \item{respLev}{a numeric vector containing the response levels.}
  \item{interval}{character string specifying the type of confidence intervals to be supplied. The default is "none". 
  Use "delta" for asymptotics-based confidence intervals (using the delta method and the t-distribution).
  Use "fls" for from logarithm scale based confidence intervals (in case the parameter in the model is log(ED50) as for
  the \code{\link{llogistic2}}) models. The only alternative for model-robust fits is using inverse regression.}
  \item{level}{numeric. The level for the confidence intervals. The default is 0.95.}
  \item{reference}{character string. Is the upper limit or the control level the reference?}
  \item{type}{character string. Whether the specified response levels are absolute or relative (default).}
  \item{lref}{numeric value specifying the lower limit to serve reference.}
  \item{uref}{numeric value specifying the upper limit to serve reference (eg. 100\%).}  
  \item{bound}{logical. If TRUE only ED values between 0 and 100\% are allowed. FALSE is useful for hormesis models.}
  \item{od}{logical. If TRUE adjustment for over-dispersion is used.}
  \item{display}{logical. If TRUE results are displayed. Otherwise they are not (useful in simulations).}
  \item{pool}{logical. If TRUE curves are pooled. Otherwise they are not. This argument only works for models with
  independently fitted curves as specified in \code{\link{drm}}.}
  \item{logBase}{numeric. The base of the logarithm in case logarithm transformed dose values are used.}  
  \item{method}{character string specifying if bisection or grid search should be used to determine ED estimates.
  Grid search may give better results for ED level close to the boundaries of the dose range or in case of
  a non-monotonous dose-response curves. The bisection method is faster than the grid search.}
  \item{cgridsize}{numeric specifying the number of grid points in the initial grid used for both bisection and grid
  search to narrow down the interval where the ED estimate is to be found.}
  \item{gridsize}{numeric specifying the number of grid points in the second finer grid search.}
  \item{lower}{numeric value specifying the lower reference limit.}
  \item{upper}{numeric specifying the upper reference limit.}  
  \item{intType}{character string specifying whether confidence or prediction intervals are requested.} 
  \item{minmax}{character string indicating if the control level should be based on the the minimum dose or the maximum
  response. The latter is more suitable for dose-response data showing hormesis.} 
  \item{n}{numeric specifying the number of simulations for the bootstrap confidence intervals.}       
  \item{seedVal}{numeric giving the seed for the random number generator used for the bootstrap confidence intervals.}   
  \item{...}{further arguments (none at the moment).}
}

\details{
  For hormesis models (\code{\link{braincousens}} and \code{\link{cedergreen}}), the additional 
  arguments \code{lower} and \code{upper} may be supplied. 
  The \code{lower} and \code{upper} arguments specify the lower and upper limits of the bisection method used to
  find the ED values.  The lower and upper limits need to be smaller/larger than the EDx level to be calculated. 
  The default limits are 0.001 and 1000 for \code{braincousens} and 0.0001 and 10000 for \code{cedergreen}, but this 
  may need to be increased. Notice that the lower limit should not be set to 0 (use instead something 
  like 1e-3, 1e-6, ...).
  
  For model-robust fits the arguments \code{lower} and \code{upper} can be used to specify reference values for the
  lower and upper limits of the dose-response relationship. This only applies for the continuous responses. For
  quantal responses, the reference values are fixed 0 and 1, respectively.
}

\value{
  A matrix with two or more columns, containing the estimates 
  and the corresponding estimated standard errors and possibly lower and upper confidence limits.
}

%\references{ ~put references to the literature/web site here ~ }
\author{Christian Ritz}

%\note{This function is only implemented for the built-in functions of class 'braincousens', 'gompertz', 'logistic' and 'mlogistic'.}

\seealso{
  The related function \code{\link{SI}}. For model-robust fits, examples are found 
  in the help of \code{\link{mrdrm}}.
}

\examples{

### How to use 'ED'

## Fitting 4-parameter log-logistic model
ryegrass.m1<-drm(ryegrass, fct = LL.4())

## Calculating EC/ED values
ED(ryegrass.m1, c(10, 50, 90)) 
## first column: the estimates of ED10, ED50 and ED90
## second column: the estimated standard errors 

### How to use the argument 'ci'

## Also displaying 95% confidence intervals
ED(ryegrass.m1, c(10,50,90), interval = "delta")

## Comparing delta method and back-transformed 
##  confidence intervals for ED values

## Fitting 4-parameter log-logistic 
##  in different parameterisation (using LL2.4)
ryegrass.m2 <- drm(ryegrass, fct = LL2.4())  

ED(ryegrass.m1, c(10,50,90), interval = "fls")
ED(ryegrass.m2, c(10,50,90), interval = "delta")


### How to use the argument 'bound'

## Fitting the Brain-Cousens model
lettuce.m1 <- drm(weight ~ conc, 
data = lettuce, fct = BC.4())

### Calculating ED[-10]

# It does not work
#ED(lettuce.m1, -10)  

## Now it does work
ED(lettuce.m1, -10, bound = FALSE)  # works
ED(lettuce.m1, -20, bound = FALSE)  # works

## It does not work for another reason: ED[-30] does not exist 
#ED(lettuce.m1, -30, bound = FALSE)  

}

\keyword{models}
\keyword{nonlinear}
