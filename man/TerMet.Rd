\name{TerMet}
\alias{TerMet}
\docType{data}
\title{The effect of terbuthylazin on growth rate}
\description{
  Test on the effect of terbuthylazin on  \it{Lemna minor}, performed on an aseptic
  culture according to the OECD-guidelines.  
}
\usage{data(TerMet)}
\format{
  A data frame with 30 observations on the following 2 variables.
  \describe{
    \item{dose}{a numeric vector}
    \item{rgr}{a numeric vector}
  }
}
\details{
  Dose is \eqn{\mu l^{-1}} and rgr is the relative growth rate of \it{Lemna}.  
}
\source{
  Cedergreen N. (2004) Unpublished bioassay data.
}
\examples{

TerMet[1:5,]  # displaying first 5 rows of the data set

}
\keyword{datasets}
