\name{braincousens}
\alias{braincousens}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{The Brain-Cousens modification of the logistic function}
\description{
  'braincousens' provides a very general way of specifying Brain-Cousens' 
  modified logistic dose response functions, under various constraints on the parameters.
}
\usage{
  braincousens(lowerc = c(-Inf, -Inf, -Inf, -Inf, -Inf), 
  upperc = c(Inf, Inf, Inf, Inf, Inf), fixed = c(NA, NA, NA, NA, NA), 
  names = c("b", "c", "d", "e", "f"), scaleDose = TRUE, useDer = FALSE)
}
\arguments{
  \item{lowerc}{numeric vector. The lower bound on parameters. Default is minus infinity.}
  \item{upperc}{numeric vector. The upper bound on parameters. Default is plus infinity.}
  \item{fixed}{numeric vector. Specifies which parameters are fixed and at what value they are fixed. NAs for parameter that are not fixed.}
  \item{names}{a vector of character strings giving the names of the parameters (should not contain ":"). The default is reasonable (see under 'Usage').
               The order of the parameters is: b, c, d, e, f (see under 'Details').}
  \item{scaleDose}{logical. If TRUE dose values are scaled around 1 during estimation; 
                   this is required for datasets where all dose values are small. }               
  \item{useDer}{logical. If TRUE derivatives are supplied, otherwise they are not supplied. Not yet implemented!}
}
\details{
  The Brain-Cousens model is given by the expression
  \deqn{ f(x) = c + \frac{d-c+fx}{1+\exp(b(\log(x)-\log(e)))}}
  which is a five-parameter model.
  
  It is a modification of the four-parameter logistic curve to take hormesis into account.
}
\value{
  The value returned by the 'braincousens' is a list with the following components
  \item{fct}{The dose response function.}
  \item{ssfct}{The self starter function.}
  \item{deriv1}{The first derivative.}  
  \item{deriv2}{The second derivative.}  
  \item{lowerc}{The lower bounds on the parameters.}
  \item{upperc}{The upper bounds on the parameters.}
  \item{edfct}{The ED function.}
  \item{sifct}{The SI function.}    
  \item{maxfct}{The function for calculating the mean maximum of the dose response curve.}  
}
\references{Brain, P. and Cousens, R. (1989) An equation to describe dose responses where there is stimulation of growth at low doses, 
            \emph{Weed Research}, \bold{29}, 93--96.}
\author{Christian Ritz}
\note{This function is for use with the function \code{\link{multdrc}}.}

\seealso{Special cases of the function 'braincousens' are \code{\link{bcl3}} and \code{\link{bcl4}}.}
\examples{

## Brain-Cousens model with the constraint b>1
model1 <- multdrc(FA, fct=braincousens(fixed=c(NA, NA, NA, NA, NA), 
lowerc=c(1, -Inf, -Inf, -Inf, -Inf)), control=mdControl(constr=TRUE))

summary(model1)

## Brain-Cousens model with the constraint f>0
model2 <- multdrc(FA, fct=braincousens(fixed=c(NA, NA, NA, NA, NA), 
lowerc=c(-Inf, -Inf, -Inf, -Inf, 0)), control=mdControl(constr=TRUE))

summary(model2)
rm(model1, model2)

}
\keyword{models}
\keyword{nonlinear}
