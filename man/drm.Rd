\name{drm}

\alias{drm}

\title{Fitting dose-response models}

\description{
  A general tools for fitting concentration/dose/time-response models.
}

\usage{
  drm(formula, curve, pmodels, weights, data = NULL, subset, fct, 
  adjust = c("none", "bc1", "bc2", "vp"), bc = NULL, bcAdd = 0,
  type = c("continuous", "binomial", "Poisson", "survival"), 
  start, start2, na.action = na.fail, hetvar = NULL, robust = "mean", logDose = NULL, 
  fctList = NULL, control = drmc(), lowerl = NULL, upperl = NULL)
}

\arguments{
  \item{formula}{a symbolic description of the model to be fit. Either of the form 'response \eqn{~} dose' 
  or as a data frame with response value in first column and dose in second column.}
  \item{curve}{a numeric vector or factor containing the grouping of the data.}
  \item{pmodels}{a data frame with a many columns as there are parameters in the non-linear function. 
  Or a list containing a formula for each parameter in the non-linear function.}
  \item{weights}{a numeric vector containing weights.}
  \item{data}{an optional data frame containing the variables in the model.}
  \item{subset}{an optional vector specifying a subset of observations to be used in the fitting process.}
  \item{fct}{a list with three or 5 elements specifying the non-linear 
  function, the accompanying self starter function, the names of the parameter in the non-linear f
  unction and, optionally, the first and second derivatives.}
  \item{adjust}{a character string specifying the type of adjustment for variance inhomogeneity.}
  \item{bc}{a numeric value specifying the lambda parameter to be used in the Box-Cox transformation.}
  \item{bcAdd}{a numeric value specifying the constant to be added on both sides prior to Box-Cox transformation. 
  The default is 0.}
  \item{type}{a character string specifying the data type: continuous is the only option currently.}
  \item{start}{an optional numeric vector containing starting values for all mean parameters in the model. 
  Overrules any self starter function.}
  \item{start2}{an optional numeric vector containing starting values for all variance parameters in the model.}
  \item{na.action}{a function which indicates what should happen when the data contain 'NA's.  The default is 'na.fail'. To omit 'NA's use
  'na.omit'.}
  \item{hetvar}{a vector specifying the grouping for heterogeneous variances.}
  \item{robust}{a character string specifying the rho function for robust estimation. Default is non-robust 
  least squares estimation ("mean"). Available robust methods are: median estimation ("median"), 
  least median of squares ("lms"), least trimmed squares ("lts"), metric trimming ("trimmed"), 
  metric winsorizing ("winsor") and Tukey's biweight ("tukey").}
  \item{logDose}{a numeric value or NULL. If log doses value are provided the base of the logarithm should be specified (exp(1) for the natural logarithm
and 10 for 10-logarithm).}
  \item{fctList}{a list of functions in case different functions need to be fit to different curves.}
  \item{control}{a list of arguments controlling constrained optimisation (zero as boundary), maximum number of iteration in the optimisation,
  relative tolerance in the optimisation, warnings issued during the optimisation.}
  \item{lowerl}{ ~~Describe \code{method} here~~ }  
  \item{upperl}{ ~~Describe \code{method} here~~ }
}

\details{
  This function is intended to supersede the function \code{\link{multdrc}}.

  This function uses the function 'optim' for the minimisation. The control arguments are specified using 
  the function \code{\link{drmc}}.
  
  If the argument \code{adjust} is "bc1" it may be necessary to constrain the optimisation 
  (using \code{\link{mdControl}}) in order to comply with the requirement of non-negative values for the 
  Box-Cox transformation. Alternatively a constant can be added in the Box-Cox transformation (using \code{bcAdd}). The residuals available will be based on the transformed response values.
  
  If \code{adjust} is "vp" then \code{robust} is set to "mean".
  
  For robust estimation MAD (median abslolute deviance) is used to estimate the residual variance.  
}

\value{
  An object of class 'drc'. 
}

%\references{ ~put references to the literature/web site here ~ }

\author{Christian Ritz and Jens C. Streibig}

\note{
  The \code{curve} argument is subsequently used for plotting dose-response curves. For large datasets 
  it need not represent the actual curves, but it could represent several curves having the parameters 
  (to avoid too many curves in the plot). The same applies for use with the function \code{\link{ED}} 
  and \code{\link{SI}}.

  The columns of a data frame argument to \code{pmodels} are automatically converted into factors. 
  This does not happen if a list is specified. The data frame argument allows collapsing parameters for 
  some or all curves (using \code{\link{colFct}}), whereas the list argument allows specification of 
  additive models for the parameters.
}

\seealso{The function \code{\link{nls}} can also be used for nonlinear regression.}

\examples{

m1 <- drm(ryegrass, fct = LL.4())
m2 <- drm(SLOPE~DOSE, CURVE, data=PestSci, fct = LL.4())
m3 <- drm(rgr~dose, pct, data=TM, fct = LL.4())

}
\keyword{models}
\keyword{nonlinear}
