\name{glymet}

\alias{glymet}

\docType{data}

\title{Glyphosate and metsulfuron-methyl tested on algae.}

\description{
  The dataset has 7 mixtures, 8 dilutions, two replicates and 5 common control controls. 
  Four observations are missing, giving a total of 113 observations.
}

\usage{data(glymet)}

\format{
  A data frame with 113 observations on the following 3 variables.
  \describe{
    \item{\code{dose}}{a numeric vector of dose values}
    \item{\code{pct}}{a numeric vector denoting the grouping according to the mixtures percentages}
    \item{\code{rgr}}{a numeric vector of response values (relative growth rates)}
  }
}

\details{
  The dataset is analysed in Soerensen et al (2007). 
  The concentration addition model can be entertained for this dataset.
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
glymet.mfree <- drm(rgr~dose, pct, data = glymet, 
fct = LL.3(), pmodels = list(~factor(pct),~1,~factor(pct) - 1),
adjust = "bc1") 


## Lack-of-fit test
anova(glymet.mfree)  # not really acceptable
summary(glymet.mfree)

## Plotting isobole structure
isobole(glymet.mfree, exchange=0.01)

## Fitting the concentration addition model
glymet.mca <- mixture(rgr~dose, pct, data = glymet, collapse = ~factor(pct), 
fct = LL.3(), model = "CA", boxcox = TRUE)

## Comparing to model with freely varying e parameter
anova(glymet.mca, glymet.mfree)  # borderline accepted

## Plotting isobole based on concentration addition
isobole(glymet.mfree, glymet.mca, exchange=0.01)  # acceptable fit

## Fitting the Hewlett model
glymet.mhew <- mixture(rgr~dose, pct, data = glymet, collapse = ~factor(pct), 
fct = LL.3(), model = "Hewlett", boxcox = TRUE)

## Comparing to model with freely varying e parameter
anova(glymet.mca, glymet.mhew)  
# borderline accepted
# the Hewlett model offers no improvement over concentration addition

## Plotting isobole based on concentration addition
isobole(glymet.mfree, glymet.mhew, exchange=0.01)  
# no improvement over concentration addition

}

\keyword{datasets}
