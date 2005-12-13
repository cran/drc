\name{g4}
\alias{g4}
\title{The four-parameter Gompertz function}
\description{
  'g4' provides the four-parameter Gompertz function, self starter function and names of the parameters.  
}
\usage{
g4(names = c("b", "c", "d", "e"), useDer = FALSE)
}
\arguments{
  \item{names}{a vector of character strings giving the names of the parameters. The default is reasonable.}
  \item{useDer}{logical. If TRUE derivatives are supplied, otherwise they are not supplied. Not implemented!}
}
\details{
  The four-parameter Gompertz model is given by the expression
  \deqn{ f(x) = c + (d-c)\exp(-\exp(b(\log(x)-e))).}
  
  The function is asymmetric about the inflection point, that is the parameter \eqn{\exp(e)}.
}
\value{
  See \code{\link{gompertz}}.
}
\references{Seber, G. A. F. and Wild, C. J (1989) \emph{Nonlinear Regression}, New York: Wiley \& Sons (pp. 330--331).}
\author{Christian Ritz}
\note{This function is for use with the function \code{\link{multdrc}}.}

\seealso{Setting \eqn{c=0} yields \code{\link{g3}}. A more flexible function, allowing fixing or constraining parameters, is \code{\link{gompertz}}.}
\examples{

model1 <- multdrc(rgr~dose, pct, data=TM, fct=g4(), cm="999")
summary(model1)

rm(model1)

}
\keyword{models}
\keyword{nonlinear}
