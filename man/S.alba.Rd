\name{S.alba}

\alias{S.alba}

\docType{data}

\title{Potency of two herbicides}

\description{
  Data are from an experiment, comparing the potency of the two herbicides glyphosate and bentazone in
  \emph{Sinapis alba}.
}

\usage{data(S.alba)}

\format{
  A data frame with 68 observations on the following 3 variables.
  \describe{
    \item{\code{Dose}}{a numeric vector giving the dose in g/ha}
    \item{\code{Herbicide}}{a factor with levels \code{Bentazone} \code{Glyphosate}}
    \item{\code{DryMatter}}{a numeric vector containing the response}
  }
}

\details{
  The lower and upper limits for the two herbicides can be assumed identical, whereas slopes and ED50 values 
  are different (in the log-logistic model).
}

\source{
  Christensen, M. G. and Teicher, H. B. and Streibig, J. C. (2003) Linking fluorescence 
  induction curve and biomass in herbicide screening, \emph{Pest Management Science},
  \bold{59},  1303--1310.

}

\seealso{
  See the examples in the examples sections in \code{\link{multdrc}} and \code{\link{SI}}.
}

\examples{

## Fitting a log-logistic model with
##  common lower and upper limits
model1 <- drm(DryMatter~Dose, Herbicide, data=S.alba, fct = LL.4(),
pmodels=data.frame(Herbicide,1,1,Herbicide), adjust="bc1") 

## Plotting fitted curve and data
plot(model1)

}
\keyword{datasets}
