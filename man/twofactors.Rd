\name{twofactors}

\alias{twofactors}

\docType{data}

\title{Dose-response curves for multiple factors}

\description{
  The dataset consists of several dose response curves for each
  combination of two factors that possibly influence the 
  dose-response relationship.
}

\usage{data(twofactors)}

\format{
  A data frame with 288 observations on the following 4 variables.
  \describe{
    \item{\code{date}}{a factor with levels \code{02-05-2005} \code{20-04-2005} \code{25-04-2005}}
    \item{\code{treatment}}{a factor with levels \code{yes} \code{no}}
    \item{\code{conc}}{a numeric vector of concentrations/doses}
    \item{\code{length}}{a numeric vector of plant lengths}
  }
}

\details{
  It turns out that the two factors only need to enter
  the model for the parameter e, in an additive way.
}

\source{
  Artificial data
}

\references{
}

\examples{

model1 <- multdrc(length ~ conc,
          date:treatment,
          data = twofactors, fct = l3(),
          collapse =
          list(~date:treatment-1, ~date:treatment-1,
               ~date:treatment-1),
          na.action=na.omit)
          
summary(model1)

rm(model1)
}

\keyword{datasets}
