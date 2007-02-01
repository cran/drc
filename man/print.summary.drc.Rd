\name{print.summary.drc}

\alias{print.summary.drc}

\title{Printing summary of non-linear model fits}

\description{
  'print.summary' displays the summary.  
}

\usage{

  \method{print}{summary.drc}(x, ...)

}

\arguments{
  \item{x}{an object of class 'drc'.}
  \item{...}{additional arguments.}
}

\value{
  'print.summary.drc' produces a formatted output of the summary statistics.
}

\author{Christian Ritz}

\examples{

## Invoking the print.summary method for objects of class 'drc'
summary(multdrc(ryegrass))  

}
\keyword{models}
\keyword{nonlinear}
