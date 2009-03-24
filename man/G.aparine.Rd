\name{G.aparine}

\alias{G.aparine}

\docType{data}

\title{Herbicide applied to Galium aparine}

\description{
  Small plants of \emph{Galium aparine}, growing in pots in a green house, were sprayed with the technical 
  grade phenmidipham herbicide either alone or in mixture with an ester of oleic acid. 
  The plants were allowed to grow in the green house for 14 days after herbicide treatment. 
  Then the dry matter was measured per pot.
}

\usage{data(G.aparine)}

\format{
  A data frame with 240 observations on the following 3 variables.
  \describe{
    \item{\code{dose}}{a numeric vector of dose value (g/ha)}
    \item{\code{drymatter}}{a numeric vector of dry matter weights (mg/pot)}
    \item{\code{treatment}}{a numeric vector giving the grouping: 0: control, 1,2: herbicide formulations}
  }
}

%\details{}

\source{
  Cabanne, F., Gaudry, J. C. and Streibig, J. C. (1999) Influence of alkyl oleates on efficacy 
  of phenmedipham applied as an acetone:water solution on Galium aparine, 
  \emph{Weed Research}, \bold{39}, 57--67.
}

%\references{}

\examples{

## Fitting a model with a common control
c.m1 <- drm(drymatter ~ dose, treatment, data = G.aparine, 
pmodels = data.frame(treatment,treatment,1,treatment), fct = LL.4())

## Visual inspection of fit
plot(c.m1, broken = TRUE)

## Lack of fit test
modelFit(c.m1)

## Summary output
summary(c.m1)

## Predicted values with se and confidence intervals
predict(c.m1, interval = "confidence")

## Calculating the relative potency
SI(c.m1, c(50,50))

## Showing the relative potency as a
## function of the response level
relpot(c.m1)
relpot(c.m1, ci = "delta")
# appears constant!

## Response level in percent
relpot(c.m1, scale="percent")

## Fitting a reduced model
c.m2 <- drm(drymatter ~ dose, treatment, data = G.aparine, 
pmodels = data.frame(1,treatment,1,treatment), fct = LL.4())

anova(c.m2, c.m1)

## Showing the relative potency 
relpot(c.m2)

## Fitting the same model in a different parameterisation
c.m3<-drm(drymatter~dose, treatment, data=G.aparine, 
pmodels=data.frame(treatment,treatment,1,treatment),fct=LL2.4())

SI(c.m3, c(50,50), logBase = exp(1))

}
\keyword{datasets}
