\name{mecter}

\alias{mecter}

\docType{data}

\title{Mechlorprop and terbythylazine tested on Lemna minor}

\description{
  Data consist of 5 mixture, 6 dilutions, three replicates, and 12 common controls; in total 102 onservations.
}

\usage{data(mecter)}

\format{
  A data frame with 102 observations on the following 3 variables.
  \describe{
    \item{\code{dose}}{a numeric vector of dose values}
    \item{\code{pct}}{a numeric vector denoting the grouping according to the mixtures percentages}
    \item{\code{rgr}}{a numeric vector of response values (relative growth rates)}
  }
}

\details{
  The dataset is analysed in Soerensen et al (2007). 
  The asymmetric Voelund model is appropriate, whereas the symmetric Hewlett model is not.
}

\source{
  The dataset is kindly provided by Nina Cedergreen, Department of Agricultural Sciences, 
  Royal Veterinary and Agricultural University, Denmark.
}

\references{
  Soerensen, H. and Cedergreen, N. and Skovgaard, I. M. and Streibig, J. C. (2007) 
  An isobole-based statistical model and test for synergism/antagonism in binary mixture toxicity experiments, 
  \emph{Environmental and Ecological Statistics}, \bold{14}, 383--397.
}

\examples{

## Fitting the model with freely varying ED50 values
mecter.mfree <- drm(rgr~dose, pct, data=mecter, 
fct=LL.4(), pmodels=list(~1,~1,~1,~factor(pct) - 1)) 

## Lack-of-fit test
anova(mecter.mfree)  # not really acceptable
summary(mecter.mfree)

## Plotting isobole structure
isobole(mecter.mfree, exchange=0.02)

## Fitting the concentration addition model
mecter.mca <- mixture(rgr~dose, pct, data=mecter, collapse=~1, 
fct=LL.4(), model="CA")

## Comparing to model with freely varying e parameter
anova(mecter.mca, mecter.mfree)  # rejected

## Plotting isobole based on concentration addition
isobole(mecter.mfree, mecter.mca, exchange=0.02)  # poor fit

## Fitting the Hewlett model
mecter.mhew <- mixture(rgr~dose, pct, data=mecter, collapse=~1, 
fct=LL.4(), model="Hewlett")

## Comparing to model with freely varying e parameter
anova(mecter.mhew, mecter.mfree)  # rejected

## Plotting isobole based on the Hewlett model
isobole(mecter.mfree, mecter.mhew, exchange=0.02)  # poor fit

## Fitting the Voelund model
mecter.mvoe<-mixture(rgr~dose, pct, data=mecter, collapse=~1, 
fct=LL.4(), model="Voelund")

## Comparing to model with freely varying e parameter
anova(mecter.mvoe, mecter.mfree)  # accepted

## Plotting isobole based on the Voelund model
isobole(mecter.mfree, mecter.mvoe, exchange=0.02)  # good fit


}
\keyword{datasets}
