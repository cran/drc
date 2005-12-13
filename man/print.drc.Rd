\name{print.drc}
\alias{print.drc}
\title{Printing key features}
\description{
  'print' displays brief information on an object of class 'drc'.
}
\usage{
print.drc(x, ..., digits = max(3, getOption("digits") - 3))
}
\arguments{
  \item{x}{an object of class 'drc'.}
  \item{...}{additional arguments.}
  \item{digits}{an integer giving the number of digits of the parameter coefficients. Default is 3.} 
}
\author{Christian Ritz}
\examples{

model1<-multdrc(FA)
print(model1)
model1  # gives the same output as the previous line

rm(model1)

}
\keyword{models}
\keyword{nonlinear}
