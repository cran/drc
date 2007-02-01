\name{CRS.5a}

\alias{CRS.5a}
\alias{CRS.5b}
\alias{CRS.5c}

\alias{ml4a}
\alias{ml4b}
\alias{ml4c}

\alias{UCRS.5a}
\alias{UCRS.5b}
\alias{UCRS.5c}

\alias{uml4a}
\alias{uml4b}
\alias{uml4c}


\title{The Cedergreen-Ritz-Streibig model}

\description{
  'CRS.5a', 'CRS.5b' and 'CRS.5c' provide the Cedergreen-Ritz-Streibig modified log-logistic model for describing
  (inverse u-shaped) hormesis.
  
  'UCRS.5a', 'UCRS.5b' and 'UCRS.5c' provide the Cedergreen-Ritz-Streibig modified log-logistic model for
  describing u-shaped hormesis. 
}

\usage{
  CRS.5a(names = c("b", "c", "d", "e", "f"))

  UCRS.5a(names = c("b", "c", "d", "e", "f"))
}

\arguments{
  \item{names}{a vector of character strings giving the names of the parameters. The default is reasonable (see above).}
}

\details{
The model is given by the expression
  \deqn{ f(x) = c + \frac{d-c+f \exp(-1/x^{\alpha})}{1+\exp(b(\log(x)-\log(e)))}}
which is a five-parameter model.
  
It is a modification of the four-parameter logistic curve to take hormesis into account.
  
The u-shaped model is given by the expression
  \deqn{ f(x) = c + d - \frac{d-c+f \exp(-1/x^{\alpha})}{1+\exp(b(\log(x)-\log(e)))}}
  
The a,b,c models are obtained by setting alpha equal to 1, 0.5 and 0.25, respectively.  
}

\value{
  See \code{\link{cedergreen}}.
}

\references{
  See the reference under \code{\link{cedergreen}}.
}

\author{Christian Ritz}

\note{
  This function is for use with the function \code{\link{drm}} or \code{\link{multdrc}}.  
}

\seealso{
  Similar functions are \code{\link{CRS.4a}} and \code{\link{UCRS.4a}}, but with the
  lower limit fixed at 0 (one parameter less).
}

\examples{

## Modified logistic model
model1 <- multdrc(hormesis[,c(2,1)], fct=CRS.5a())
summary(model1)
ED(model1, c(50))

model2 <- multdrc(hormesis[,c(2,1)], fct=CRS.5b())
summary(model2)
ED(model2, c(50))

model3 <- multdrc(hormesis[,c(2,1)], fct=CRS.5c())
summary(model3)
ED(model3, c(50))

rm(model1, model2, model3)

}
\keyword{models}
\keyword{nonlinear}
