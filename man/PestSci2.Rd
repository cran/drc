\name{PestSci2}
\alias{PestSci2}
\docType{data}
\title{Bioassay data with two factors}
\description{
  The data set consists of 4 dose response curves. For each dose and response, the level of two factors was recorded.
}
\usage{data(PestSci2)}
\format{
  A data frame with 84 observations on the following 5 variables.
  \describe{
    \item{CURVE}{a numeric vector}
    \item{DOSE}{a numeric vector}
    \item{SLOPE}{a numeric vector}
    \item{B}{a factor with levels \code{1} \code{2}}
    \item{A}{a factor with levels \code{a} \code{b}}
  }
}
\details{
  A constructed data set.
}
\examples{

PestSci2[1:10,]  # displaying first 10 rows

}
\keyword{datasets}
