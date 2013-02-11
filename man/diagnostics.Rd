\name{diagnostics}

\alias{diagnostics}

\title{Information on estimation}

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
  The function displays the convergence status and the number of evaluations used.
}

\author{Christian Ritz}

\examples{

terbuthylazin.m1 <- drm(rgr~dose, data = terbuthylazin, fct = LL.4())
diagnostics(terbuthylazin.m1)

terbuthylazin.m2 <- drm(rgr~dose, data = terbuthylazin, start=c(1,0.01,0.38,130), fct = LL.4()) 
diagnostics(terbuthylazin.m2)  
## the start values given result in faster convergence

}
\keyword{models}
\keyword{nonlinear}
