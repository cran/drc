\name{finney71}

\alias{finney71}

\docType{data}

\title{Example from Finney (1971)}

\description{
  For each of six concentration of an insecticid the number of insects affected (out of the number of insects)
  was recorded.
}

\usage{data(finney71)}

\format{
  A data frame with 6 observations on the following 3 variables.
  \describe{
    \item{\code{dose}}{a numeric vector}
    \item{\code{total}}{a numeric vector}
    \item{\code{affected}}{a numeric vector}
  }
}

%\details{}

\source{
   Finney, D. J. (1971) \emph{Probit Analysis}, Cambridge: Cambridge University Press.
}

%\references{}

\examples{

## Fitting a logistic regression model
model1 <- multdrc(affected/total~dose, weights=total, data=finney71,
fct=l2(), type="binomial")
summary(model1)
plot(model1)


rm(model1)
}
\keyword{datasets}
