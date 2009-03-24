\name{S.alba}

\alias{S.alba}

\docType{data}

\title{Potency of two herbicides}

\description{
  Data are from an experiment, comparing the potency of the two herbicides glyphosate and bentazone in
  white mustard \emph{Sinapis alba}.
}

\usage{data(S.alba)}

\format{
  A data frame with 68 observations on the following 3 variables.
  \describe{
    \item{\code{Dose}}{a numeric vector containing the dose in g/ha.}
    \item{\code{Herbicide}}{a factor with levels \code{Bentazone} \code{Glyphosate} (the two herbicides applied).}
    \item{\code{DryMatter}}{a numeric vector containing the response (dry matter in g/pot).}
  }
}

\details{
  The lower and upper limits for the two herbicides can be assumed identical, whereas slopes and ED50 values 
  are different (in the log-logistic model).
}

\source{
  Christensen, M. G. and Teicher, H. B., and Streibig, J. C. (2003) Linking fluorescence 
  induction curve and biomass in herbicide screening, \emph{Pest Management Science},
  \bold{59},  1303--1310.

}

\seealso{
  See the examples in the examples sections in \code{\link{drm}} and \code{\link{SI}}.
}

\examples{

## Fitting a log-logistic model with
##  common lower and upper limits
S.alba.m1 <- drm(DryMatter~Dose, Herbicide, data=S.alba, fct = LL.4(),
pmodels=data.frame(Herbicide,1,1,Herbicide)) 

## Applying the optimal transform-both-sides Box-Cox transformation
S.alba.m2 <- boxcox(S.alba.m1, method = "anova") 

## Plotting fitted regression curves together with the data
plot(S.alba.m2)

}
\keyword{datasets}
