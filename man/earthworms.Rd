\name{earthworms}

\alias{earthworms}

\docType{data}

\title{Earthworms}

\description{
  Dataset obtained from a toxicity test using earthworms, counting the number of earthworms staying in the toxic earth container (not migrating to the neighbouring container of clean earth).
}

\usage{data(earthworms)}

\format{
  A data frame with 35 observations on the following 3 variables.
  \describe{
    \item{\code{dose}}{a numeric vector}
    \item{\code{number}}{a numeric vector}
    \item{\code{total}}{a numeric vector}
  }
}

\details{
  At dose 0 around half of the earthworms wil be in each of the two containers. Thus it is not appropriate to this an ordinary logistic regression (with upper limit 1). 
}

\source{
  The dataset is kindly provided by Nina Cedergreen, Royal Veterinary and Agricultural University, Denmark.
}

%\references{}

\examples{

## Fitting a logistic regression model
earthworms.m1 <- drm(number/total~dose, weights = total, data = earthworms,
fct = LL.2(), type = "binomial")
modelFit(earthworms.m1)  # a crude goodness-of-fit test

## Fitting an extended logistic regression model 
##  where the upper limit is estimated
earthworms.m2 <- drm(number/total~dose, weights = total, data = earthworms,
fct = LL.3(), type = "binomial")
modelFit(earthworms.m2)  # goodness-of-fit test
# improvement not visible in test!!!

## Comparing model1 and model2 (can model2 be reduced to model1?)
anova(earthworms.m1, earthworms.m2, test = "Chisq")

}
\keyword{datasets}
