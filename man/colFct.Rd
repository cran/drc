\name{colFct}
\alias{colFct}
\title{Collapse of levels}
\description{
  Collapse levels in a vector into fewer levels. 
}
\usage{
colFct(assayNo, noVec)
}
\arguments{
  \item{assayNo}{a vector to be collapsed.}
  \item{noVec}{a vector of indices in \code{assayNo} that should be collapsed.}
}
\details{
  This function is primarily for use in the collapse argument in \code{\link{multdrc}}.
}
\value{
  A vector of same length as \code{assayNo} containing the new levels.
}
\author{Christian Ritz}
%\note{ ~~further notes~~ }

%\seealso{See the example section for \code{\link{multdrc}}.}
\examples{

x <- c("a","a","b","b","b","b","c","c","d","d","d","e")
colFct(x,c(4,5))
rm(x)

}
\keyword{models}
\keyword{nonlinear}
