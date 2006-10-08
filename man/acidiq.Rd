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
  A data frame with 180 observations on the following 4 variables.
  \describe{
    \item{\code{cell}}{a numeric vector}
    \item{\code{dose1}}{a numeric vector}
    \item{\code{pct1}}{a numeric vector}
    \item{\code{rgr1}}{a numeric vector}
  }
}

\details{
  The dataset is analysed in S{\o}rensen et al (2006). Hewlett's symmetric model seems appropriate for this dataset. 
}

\source{
  The dataset is kindly provided by Nina Cedergreen, Department of Agricultural Sciences, Royal Veterinary and Agricultural University, Denmark.
}

\references{
  S{\o}rensen, H. and Cedergreen, N. and Skovgaard, I. M. and Streibig, J. C. (2006) An isobole-based statistical model and test 
  for synergism/antagonism in binary mixture toxicity experiments, \emph{Statistical Ecology and Environmental Statistics}, \bold{13}. ?--??.
}

\examples{

## Fitting the model with freely varying ED50 values
model1 <- multdrc(rgr1~dose1, pct1, data=acidiq, 
collapse=list(~factor(pct1),~1,~1,~factor(pct1) - 1), 
boxcox=TRUE)
summary(model1)


## Fitting the concentration addition model
model1ca <- mixture(rgr1~dose1, pct1, data=acidiq, collapse=~factor(pct1), 
model="CA", boxcox=TRUE)
anova(model1ca)  # comparison to model with freely varying e parameter
summary(model1ca)


## Fitting the Hewlett model
model1h <- mixture(rgr1~dose1, pct1, data=acidiq, collapse=~factor(pct1), 
model="Hewlett", boxcox=TRUE)
anova(model1h)  # comparison to a model with freely varying e parameter
anova(model1ca, model1h)  # rejecting concentration addition
summary(model1h)


## Fitting the Voelund model
model1v <- mixture(rgr1~dose1, pct1, data=acidiq, collapse=~factor(pct1), 
model="Voelund", boxcox=TRUE)
anova(model1v)  # comparison to a model with freely varying e parameter
summary(model1v)


rm(model1, model1ca, model1h, model1v)
}
\keyword{datasets}
