\name{MAX}

\alias{MAX}

\title{Maximum mean response}

\description{
  'MAX' estimates the maximum mean response and the dose at which it occurs.
}

\usage{
MAX(obj, upper=NULL)
}

\arguments{
  \item{obj}{an object of class 'drc'.}
  \item{upper}{numeric. Upper limit for bisection method. Need to be larger than EDx level to be calculated.}
}

\details{
  This function is only implemented for the built-in functions of class 'braincousens' and 'mlogistic'.
}

\value{
  A matrix with one row per curve in the data set and two columns: 
  one containing the dose at which the maximum occurs 
  and one containing the corresponding maximum response.
}

%\references{ ~put references to the literature/web site here ~ }

\author{Christian Ritz}

%\note{This function is only implemented for the built-in functions of class 'braincousens' and 'mlogistic'.}

%\seealso{The related function \code{\link{SI}}.}

\examples{

model1 <- multdrc(hormesis[,c(2,1)], fct=ml3c())
MAX(model1)

rm(model1)

}
\keyword{models}
\keyword{nonlinear}
