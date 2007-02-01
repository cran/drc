\name{ED}
\alias{ED}

\title{Estimating effective doses}

\description{
  'ED' estimates effective doses given (ECp/EDp/ICp) for given reponse levels.
}

\usage{
  ED(object, respLev, bound = TRUE, od = FALSE, ci = c("none", "delta", "fls"), 
  level = ifelse(!(ci=="none"), 0.95, NULL), logBase = NULL,
  reference = c("upper", "control"), type = c("relative", "absolute"), ...)
}

\arguments{
  \item{object}{an object of class 'drc'.}
  \item{respLev}{a numeric vector containing the response levels.}
  \item{bound}{logical. If TRUE only ED values between 0 and 100\% are allowed. FALSE is useful for hormesis models.}
  \item{od}{logical. If TRUE adjustment for over-dispersion is used.}
  \item{ci}{character string specifying the type of confidence intervals to be supplied. The default is "none". 
  Use "delta" for asymptotics-based confidence intervals (using the delta method and the t-distribution).
  Use "fls" for from logarithm scale based confidence intervals (in case the parameter in the model is log(ED50) as for
  the \code{\link{boltzmann}}) models.}
  \item{level}{numeric. The level for the confidence intervals. The default is 0.95.}
  \item{logBase}{numeric. The base of the logarithm in case logarithm transformed dose values are used.}
  \item{reference}{character string. Is the upper limit or the control level the reference?}
  \item{type}{character string. Whether the specified response levels are absolute or relative (default).}
  \item{...}{additional arguments for the function doing the calculations.}
}

\details{
  This function is only implemented for the built-in functions of class 'braincousens', 'gompertz', 'logistic' and 'mlogistic'.

  For objects of class 'braincousens' or 'mlogistic' the additional argument may be the 'upper' argument
  or the 'interval' argument. 
  The 'upper' argument specifies the upper limit of the bisection method. 
  The upper limit needs to be larger than the EDx level to be calculated. The default limit is 1000, but
  this may need to be increased.
  The 'interval' argument should specify a rough interval in which the dose
  yielding the maximum hormetical response lies. The default interval is 'c(0.001, 1000)'. Notice that the lower limit
  should not be set to 0 (use something like 1e-3, 1e-6, ...).
}

\value{
  A matrix with two columns: one containing the estimates and one containing the corresponding estimated standard errors.
}

%\references{ ~put references to the literature/web site here ~ }
\author{Christian Ritz}

%\note{This function is only implemented for the built-in functions of class 'braincousens', 'gompertz', 'logistic' and 'mlogistic'.}

\seealso{The related function \code{\link{SI}}.}

\examples{

model1<-multdrc(ryegrass)

ED(model1, c(10,50,90)) 
## first column: the estimates of ED10, ED50 and ED90
## second column: the estimated standard errors 

## Also displaying 95% confidence intervals
ED(model1, c(10,50,90), ci = "delta")

## The argument 'bound'
model2 <- multdrc(weight ~ conc, 
data = hormesis, fct = bcl3())

## Calculating ED[-10]
#ED(model2, -10)  
# does not work

ED(model2, -10, bound = FALSE)  # works
ED(model2, -20, bound = FALSE)  # works

#ED(model2, -30, bound = FALSE)  
# does not work: ED[-30] does not exist 

}

\keyword{models}
\keyword{nonlinear}
