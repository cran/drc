\name{diagnostics}
\alias{diagnostics}
\title{Display information on estimation procedure}
\description{
  Displays information on the convergence of the estimation procedure.
}
\usage{
diagnostics(object)
}
\arguments{
  \item{object}{an object of class 'drc'.}
}
\details{
  Currently this function displays the convergence status and the number of evaluations used.
}
\author{Christian Ritz}
\examples{

model1 <- multdrc(rgr~dose, data=TerMet)
diagnostics(model1)

model2 <- multdrc(rgr~dose, data=TerMet, startVal=c(1,0.01,0.38,130)) 
diagnostics(model2)  
## the start values given yield swifter convergence

rm(model1, model2)

}
\keyword{models}
\keyword{nonlinear}
