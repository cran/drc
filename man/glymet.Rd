\name{glymet}

\alias{glymet}

\docType{data}

\title{Glyphosate and metsulfuron-methyl tested on algae.}

\description{
  The dataset has 7 mixtures, 8 dilutions, two replicates and 5 common control controls. Four observations are missing, giving a total
  of 113 observations.
}

\usage{data(glymet)}

\format{
  A data frame with 113 observations on the following 4 variables.
  \describe{
    \item{\code{cell}}{a numeric vector}
    \item{\code{dose1}}{a numeric vector}
    \item{\code{pct1}}{a numeric vector}
    \item{\code{rgr1}}{a numeric vector}
  }
}

\details{
  The dataset is analysed in S{\o}rensen et al (2006). The hypothesis of concentration addition can be entertained for this dataset.
}

\source{
  The dataset is kindly provided by Nina Cedergreen, Department of Agricultural Sciences, Royal Veterinary and Agricultural University, Denmark.
}

\references{
  S{\o}rensen, H. and Cedergreen, N. and Skovgaard, I. M. and Streibig, J. C. (2006) An isobole-based statistical model and test 
  for synergism/antagonism in binary mixture toxicity experiments, \emph{Statistical Ecology and Environmental Statistics}, \bold{13}. ?--??.
}

\examples{

## Fitting the concentration addition model
model1ca <- mixture(rgr1~dose1, pct1, data=glymet, collapse=~factor(pct1), 
fct=l3(), model="CA", boxcox=TRUE)
anova(model1ca)  # comparison to model with freely varying e parameter
summary(model1ca)


## Fitting the Hewlett model
model1h <- mixture(rgr1~dose1, pct1, data=glymet, collapse=~factor(pct1), 
fct=l3(), model="Hewlett", boxcox=TRUE)
anova(model1h)  # comparison to a model with freely varying e parameter
anova(model1ca, model1h)  # concentration addition is accepted
summary(model1h)


## Fitting the Voelund model
model1v <- mixture(rgr1~dose1, pct1, data=glymet, collapse=~factor(pct1), 
fct=l3(), model="Voelund", boxcox=TRUE)
anova(model1v)  # comparison to a model with freely varying e parameter
summary(model1v)


rm(model1ca, model1h, model1v)
}

\keyword{datasets}
