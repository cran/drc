\name{bcl3}
\alias{bcl3}
\title{The three-parameter logistic function with Brain-Cousens' modification}
\description{
  'bcl3' provides the four-parameter function, self starter function and names of the parameters.
}
\usage{
bcl3(names = c("b", "d", "e", "f"), useDer = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{names}{a vector of character strings giving the names of the parameters. The default is reasonable.}
  \item{useDer}{logical. If TRUE derivatives are supplied, otherwise they are not supplied. Not implemented!}
}
\details{
  The Brain-Cousens model is given by the expression
  \deqn{ f(x) = 0 + \frac{d-0+fx}{1+\exp(b(\log(x)-\log(e)))}}
  which is a four-parameter model.
  
  It is a modification of the three-parameter logistic curve to take hormesis into account.
}
\value{
  See \code{\link{braincousens}}.
}
\references{van Ewijk, P. H. and Hoekstra, J. A. (1993) Calculation of the EC50 and its Confidence Interval When Subtoxic Stimulus Is Present,
\emph{ECOTOXICOLOGY AND ENVIRONMENTAL SAFETY}, \bold{25}, 25--32.}
\author{Christian Ritz}
\note{This function is for use with the function \code{\link{multdrc}}.}
\examples{

model1 <- multdrc(hormesis[,c(2,1)], fct=bcl3())
summary(model1)
ED(model1, c(50))  # compare the parameter estimate and 
                   # its estimated standard error 
                   # to the values in the above paper
rm(model1)

}
\keyword{models}
\keyword{nonlinear}
