\name{hormesis}

\alias{hormesis}

\docType{data}

\title{Hormesis data set}

\description{
  Results from an experiment with isobutylalcohol in nutrient solution.
}

\usage{data(hormesis)}

\format{
  A data frame with 21 observations on the following 2 variables.
  \describe{
    \item{conc}{a numeric vector}
    \item{weight}{a numeric vector}
  }
}

\details{
  The data set illustrates hormesis, presence of a subtoxic stimulus.
}

\source{
  van Ewijk, P. H. and Hoekstra, J. A. (1993) Calculation of the EC50 and its Confidence Interval When Subtoxic Stimulus Is Present,
  \emph{ECOTOXICOLOGY AND ENVIRONMENTAL SAFETY}, \bold{25}, 25--32.
}

\references{
  van Ewijk, P. H. and Hoekstra, J. A. (1994) Curvature Measures and Confidence Intervals for the Linear Logistic Model,
  \emph{Appl. Statist.}, \bold{43}, 477--487.
}

\examples{

## Look at data
head(hormesis)

## Monotonous dose-response model
hm.m1 <- drm(weight~conc, data=hormesis, fct=LL.3())

plot(hm.m1)

## Model fit in van Ewijk and Hoekstra (1994)
hm.m2 <- drm(weight~conc, data=hormesis, fct=BC.4())

plot(hm.m2, type = "fit", add = TRUE, lty = 2)

## Hormesis effect only slightly significant
summary(hm.m2)

## Hormesis effect highly significant
anova(hm.m1,hm.m2)

}
\keyword{datasets}
