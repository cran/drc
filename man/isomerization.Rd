\name{isomerization}

\alias{isomerization}

\docType{data}

\title{Reaction rate of the catalytic isomerization}

\description{
  Reaction rate of n-pentane to iso-pentane depends on various factors, 
  such as partial pressures of the products involved in the reaction.
}

\usage{data(isomerization)}

\format{
  A data frame with 24 observations on the following 4 variables.
  \describe{
    \item{\code{h}}{a numeric vector of partial pressure of hydrogen}
    \item{\code{p}}{a numeric vector of partial pressure of n-pentane}
    \item{\code{i}}{a numeric vector of partial pressures of iso-pentane}
    \item{\code{rate}}{a numeric vector of reaction rates}
  }
}

\details{
  Data are listed and thoroughly analysed in Huet et al (2004), but they have been used
  in other books as well.
}

\source{
  Huet, S., Bouvier, A., Poursat, M.-A. and Jolivet, E. (2004) 
  \emph{Statistical tools for nonlinear regression: a practical guide to S-PLUS and R}, 
  Second Edition, New York: Springer-Verlag (pp. 9--10).
}

\references{
}

\examples{

data(isomerization)

}

\keyword{datasets}
