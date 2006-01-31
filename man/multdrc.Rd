\name{multdrc}
\alias{multdrc}

\title{Fitting multiple dose-response curves}

\description{
  'multdrc' fits non-linear regression models to multiple dose-response curves.
}

\usage{
multdrc(formula, curve, collapse, weights, data=NULL, boxcox = FALSE, 
bcAdd = 0, varPower = FALSE, startVal, fct = l4(), na.action = na.fail, 
hetvar = NULL, robust="mean", type = "continuous", cm = NULL, 
logDose = NULL, control = mdControl())
}

\arguments{
  \item{formula}{a symbolic description of the model to be fit. Either of the form 'response \eqn{~} dose' or as a data frame with response value
                 in first column and dose in second column.}
  \item{curve}{a numeric vector or factor containing the grouping of the data.}
  \item{collapse}{a data frame with a many columns as there are parameters in the non-linear function. Or a list containing a formula for each 
  parameter in the non-linear function.}
  \item{weights}{a numeric vector containing weights.}
  \item{data}{an optional data frame containing the variables in the model.}
  \item{boxcox}{logical or numeric. If TRUE the optimal Box-Cox transformation is applied to the model. 
  If FALSE (the default) no Box-Cox transformation is applied. If numeric the specified value is used in the Box-Cox transformation.}
  \item{bcAdd}{a numeric value specifying the constant to be added on both sides prior to Box-Cox transformation. The default is 0.}
  \item{varPower}{logical. If TRUE the variance is modelled as a power function of the mean. 
  If FALSE (the default) no variance modelling is applied.}
  \item{startVal}{an optional numeric vector containing start values for all parameters in the model. Overrules any self starter function.}
  \item{fct}{a list with three or 5 elements specifying the non-linear 
  function, the accompanying self starter function, the names of the parameter in the non-linear function and, optionally, the first and second derivatives.}
  \item{na.action}{a function which indicates what should happen when the data contain 'NA's.  The default is 'na.fail'. To omit 'NA's use
  'na.omit'.}
  \item{hetvar}{a vector specifying the grouping for heterogeneous variances.}
  \item{robust}{a character string specifying the rho function for robust estimation. Default is non-robust least squares estimation ("mean").
  Available robust methods are: median estimation ("median"), least median of squares ("lms"), least trimmed squares ("lts"), 
  metric trimming ("trimmed"), metric winsorizing ("winsor") and Tukey's biweight ("tukey").}
  \item{type}{a character string specifying the data type: binomial and continuous are the only options currently.}
  \item{cm}{character string or numeric value specifying the level in \code{curve} corresponding to control measurements.}
  \item{logDose}{a numeric value or NULL. If log doses value are provided the base of the logarithm should be specified (exp(1) for the natural logarithm
and 10 for 10-logarithm).}
  \item{control}{a list of arguments controlling constrained optimisation (zero as boundary), maximum number of iteration in the optimisation,
  relative tolerance in the optimisation, warnings issued during the optimisation.}
}

\details{
  This function uses the function 'optim' for the minimisation. The control arguments are specified using the function
  \code{\link{mdControl}}.
  
  If the argument \code{boxcox} is TRUE it may be necessary to constrain the optimisation (using \code{\link{mdControl}}) 
  in order to comply with the requirement of non-negative values for the Box-Cox transformation. Alternatively a constant can be
  added in the Box-Cox transformation (using \code{bcAdd}). The residuals available will be based on the transformed response values.
  
  If \code{varPower} is TRUE then 'robust="mean"'. Also if "boxcox=TRUE" then 'varPower=FALSE'.
  
  For robust estimation MAD (median abslolute deviance) is used to estimate the residual variance.  
}

\value{
  An object of class 'drc'. 
}

\author{Christian Ritz and Jens C. Streibig}

\note{The \code{curve} argument is subsequently used for plotting dose-response curves. For large datasets it need not represent the actual
      curves, but it could represent several curves having the parameters (to avoid too many curves in the plot). The same applies for use
      with the function \code{\link{ED}} and \code{\link{SI}}.

      The columns of a data frame argument to \code{collapse} are automatically converted into factors. This does not happen if a list is specified.
      The data frame argument allows collapsing parameters for some or all curves (using \code{\link{colFct}}), whereas the list argument 
      allows specification of additive models for the parameters.}

\seealso{The function \code{\link{nls}} can also be used for nonlinear regression.}

\examples{

## Fitting response MEANLR to dose MM in the data set FA 
##  (type ?FA to learn about the data set)
model1 <- multdrc(MEANLR~MM, data=FA)
summary(model1)

## Fitting a 5-parameter logistic model
model2 <- multdrc(FA, fct = l5())
summary(model2)

## Fitting a 4-parameter logistic model
##  with Box-Cox transformation
model3 <- multdrc(FA, boxcox = TRUE)  
summary(model3)

## Example using 'varPower' argument 
##  modelling variance as power of mean
model4 <- multdrc(FA, varPower = TRUE)
summary(model4)

## Example using 'logDose' argument
logdose <- c(-8.045757, -7.568636, -7.086186, -6.607303, -6.130768, -5.653647, -5.180456, -4.698970)
response <- c(7.941, 16.619, 11.021, 19.767, 30.487, 64.434, 74.422, 52.472)
model5 <- multdrc(response~logdose, logDose=10)
plot(model5, conLevel=-9)  # the default level is too high

rm(model1, model2, model3, model4, model5)

}
\keyword{models}
\keyword{nonlinear}
