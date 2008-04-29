\name{acidiq}

\alias{acidiq}

\docType{data}

\title{Acifluorfen and diquat tested on Lemna minor.}

\description{
  Data from an experiment where the chemicals acifluorfen and diquat tested on Lemna minor. The dataset has 7 mixtures used in
  8 dilutions with three replicates and 12 common controls, in total 180 observations.
}

\usage{data(acidiq)}

\format{
  A data frame with 180 observations on the following 3 variables.
  \describe{
    \item{\code{dose}}{a numeric vector of dose values}
    \item{\code{pct}}{a numeric vector denoting the grouping according to the mixtures percentages}
    \item{\code{rgr}}{a numeric vector of response values (relative growth rates)}
  }
}

\details{
  The dataset is analysed in Soerensen et al (2007). 
  Hewlett's symmetric model seems appropriate for this dataset. 
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
acidiq.mfree <- drm(rgr~dose, pct, data=acidiq, fct = LL.4(),
pmodels = list(~factor(pct),~1,~1,~factor(pct) - 1), 
adjust = "bc1")

## Lack-of-fit test
anova(acidiq.mfree)
summary(acidiq.mfree)

## Plotting isobole structure
isobole(acidiq.mfree, exchange=0.1, xlim=c(0, 400), ylim=c(0,45))

## Fitting the concentration addition model
acidiq.mca <- mixture(rgr~dose, pct, data=acidiq, fct = LL.4(), 
collapse=~factor(pct), model="CA", boxcox=TRUE)

## Comparing to model with freely varying e parameter
anova(acidiq.mca, acidiq.mfree)  # rejected

## Plotting isobole based on concentration addition
isobole(acidiq.mfree, acidiq.mca, exchange=0.1, xlim=c(0, 420), ylim=c(0,45))  # poor fit

## Fitting the Hewlett model
acidiq.mhew <- mixture(rgr~dose, pct, data=acidiq, fct = LL.4(), 
collapse=~factor(pct), model="Hewlett", boxcox=TRUE)

## Comparing to model with freely varying e parameter
anova(acidiq.mfree, acidiq.mhew)  # accepted
summary(acidiq.mhew)

## Plotting isobole based on concentration addition
isobole(acidiq.mfree, acidiq.mhew, exchange=0.1, xlim=c(0, 400), ylim=c(0,45))  # good fit

}
\keyword{datasets}
