\name{PR}

\alias{PR}

\title{Expected response}

\description{
  'PR' returns the expected response for specified dose values.
}

\usage{
PR(object, xVec, ...)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{object}{ ~~Describe \code{object} here~~ }
  \item{xVec}{ ~~Describe \code{xVec} here~~ }
  \item{\dots}{ ~~Describe \code{\dots} here~~ }
}

\details{
  ~~ If necessary, more details than the description above ~~
}

\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
}

\references{ ~put references to the literature/web site here ~ }

\author{Christian Ritz after a suggestion from Andrew Kniss}

\note{ ~~further notes~~ 

 ~Make other sections like Warning with \section{Warning }{....} ~
}

\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }

\examples{

m1 <- multdrc(ryegrass)
PR(m1, c(5, 10))

m2 <- multdrc(SLOPE~DOSE, CURVE, data=PestSci)
PR(m2, c(5, 10))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{nonlinear}
