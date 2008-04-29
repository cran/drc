\name{braincousens}

\alias{braincousens}

\title{The Brain-Cousens hormesis models}

\description{
  'braincousens' provides a very general way of specifying Brain-Cousens' 
  modified log- logistic model for describing hormesis, under various constraints on the parameters.
}

\usage{
  braincousens(lowerc = c(-Inf, -Inf, -Inf, -Inf, -Inf), 
  upperc = c(Inf, Inf, Inf, Inf, Inf), fixed = c(NA, NA, NA, NA, NA), 
  names = c("b", "c", "d", "e", "f"), scaleDose = TRUE,
  fctName, fctText)
}

\arguments{
  \item{lowerc}{numeric vector. The lower bound on parameters. Default is minus infinity.}
  \item{upperc}{numeric vector. The upper bound on parameters. Default is plus infinity.}
  \item{fixed}{numeric vector. Specifies which parameters are fixed and at what value they are fixed. NAs for parameter that are not fixed.}
  \item{names}{a vector of character strings giving the names of the parameters (should not contain ":"). The default is reasonable (see under 'Usage').
               The order of the parameters is: b, c, d, e, f (see under 'Details').}
  \item{scaleDose}{logical. If TRUE dose values are scaled around 1 during estimation; 
                   this is required for datasets where all dose values are small. } 
  \item{fctName}{optional character string used internally by convenience functions.}  
  \item{fctText}{optional character string used internally by convenience functions.}                            
}

\details{
  The Brain-Cousens model is given by the expression
  \deqn{ f(x) = c + \frac{d-c+fx}{1+\exp(b(\log(x)-\log(e)))}}
  which is a five-parameter model.
  
  It is a modification of the four-parameter logistic curve to take hormesis into account.
}

\value{
  The value returned is a list containing the non-linear function, the self starter function, 
  the parameter names and additional model specific objects.  
}

\references{
  Brain, P. and Cousens, R. (1989) An equation to describe dose responses 
  where there is stimulation of growth at low doses, 
  \emph{Weed Research}, \bold{29}, 93--96.
}
            
\author{Christian Ritz}

\note{This function is for use with the function \code{\link{drm}} or \code{\link{multdrc}}.}

\seealso{The convenience functions of \code{braincousens} are \code{\link{BC.4}} and \code{\link{BC.5}}.}

\examples{

## Brain-Cousens model with the constraint b>1
model1 <- multdrc(ryegrass, fct=braincousens(fixed=c(NA, NA, NA, NA, NA), 
lowerc=c(1, -Inf, -Inf, -Inf, -Inf)), control=mdControl(constr=TRUE))

summary(model1)

## Brain-Cousens model with the constraint f>0
model2 <- multdrc(ryegrass, fct=braincousens(fixed=c(NA, NA, NA, NA, NA), 
lowerc=c(-Inf, -Inf, -Inf, -Inf, 0)), control=mdControl(constr=TRUE))

summary(model2)
rm(model1, model2)

}

\keyword{models}
\keyword{nonlinear}

\concept{hormesis hormetic effect initial stimulation}
