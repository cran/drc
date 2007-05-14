\name{mixture}

\alias{mixture}

\title{Fitting mixture models}

\description{
  'mixture' fits a concentration addition, Hewlett or Voelund model to data from mixture toxicity experiments.
}

\usage{
mixture(formula, curve, collapse, weights, data = NULL, boxcox = FALSE, 
bcAdd = 0, varPower = FALSE, startVal, fct = l4(), na.action = na.fail, 
robust = "mean", type = "continuous", cm = NULL, logDose = NULL, 
control = mdControl(), model = "Hewlett", startVal2)
}

\arguments{
  \item{formula}{a symbolic description of the model to be fit. Either of the form 'response \eqn{~} dose' or as a data frame with response value
                 in first column and dose in second column.}
  \item{curve}{a numeric vector (percentages!) containing the grouping of the data.}
  \item{collapse}{a formula starting with a tilde, specifying the model for the b parameter in the logistic model (the only 'free' parameter).}
  \item{weights}{a numeric vector containing weights.}
  \item{data}{a data frame containing the variables in the model. This argument is not optional!}
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
  \item{robust}{a character string specifying the rho function for robust estimation. Default is non-robust least squares estimation ("mean").
  Available robust methods are: median estimation ("median"), least median of squares ("lms"), least trimmed squares ("lts"), 
  metric trimming ("trimmed"), metric winsorizing ("winsor") and Tukey's biweight ("tukey").}
  \item{type}{a character string specifying the data type: binomial and continuous are the only options currently.}
  \item{cm}{character string or numeric value specifying the level in \code{curve} corresponding to control measurements.}
  \item{logDose}{a numeric value or NULL. If log doses value are provided the base of the logarithm should be specified (exp(1) for the natural logarithm
and 10 for 10-logarithm).}
  \item{control}{a list of arguments controlling constrained optimisation (zero as boundary), maximum number of iteration in the optimisation,
  relative tolerance in the optimisation, warnings issued during the optimisation.}  
  \item{model}{a character string. It can be "CA", "Hewlett" or "Voelund".}
  \item{startVal2}{an optional numeric vector supplying the lambda parameter in the Hewlett model or the eta parameters (two parameters) 
  in the Voelund model.}
}

\details{
  The function is a wrapper to \code{\link{multdrc}}, implementing the models described in S{\o}rensen et al. (2007). 
  See the paper for a discussion of the merits of the different models.
}

\value{
  An object of class 'drc'. 
}

\references{
  S{\o}rensen, H. and Cedergreen, N. and Skovgaard, I. M. and Streibig, J. C. (2007) 
  An isobole-based statistical model and test for synergism/antagonism in binary mixture toxicity experiments, 
  \emph{Statistical Ecology and Environmental Statistics}, \bold{14}. ???--???.
}

\author{Christian Ritz}

\seealso{The function has argument similar to those in \code{\link{multdrc}}.}

\examples{

## See the example under ?acidiq

}
\keyword{models}
\keyword{nonlinear}
