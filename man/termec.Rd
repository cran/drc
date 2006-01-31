\name{termec}

\alias{termec}

\docType{data}

\title{Mechlorprop and terbythylazine tested on Lemna minor}

\description{
  Data consist of 5 mixture, 6 dilutions, three replicates, and 12 common controls; in total 102 onservations.
}

\usage{data(termec)}

\format{
  A data frame with 102 observations on the following 4 variables.
  \describe{
    \item{\code{cell}}{a numeric vector}
    \item{\code{dose1}}{a numeric vector}
    \item{\code{pct1}}{a numeric vector}
    \item{\code{rgr1}}{a numeric vector}
  }
}

\details{
  The dataset is analysed in S{\o}rensen et al (2006). The asymmetric Voelund model is appropriate, whereas the symmetric Hewlett model is not.
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
model1ca <- mixture(rgr1~dose1, pct1, data=termec, collapse=~1, 
fct=l3(), model="CA", boxcox=TRUE, bcAdd=0.02)
anova(model1ca)  # comparison to model with freely varying e parameter
summary(model1ca)


## Fitting the Hewlett model
model1h <- mixture(rgr1~dose1, pct1, data=termec, collapse=~1, 
fct=l3(), model="Hewlett", boxcox=TRUE, bcAdd=0.02)
anova(model1h)  # comparison to a model with freely varying e parameter
anova(model1ca, model1h)
summary(model1h)


## Fitting the Voelund model
model1v<-mixture(rgr1~dose1, pct1, data=termec, collapse=~1, 
fct=l3(), model="Voelund", boxcox=TRUE, bcAdd=0.02)
anova(model1v)  # comparison to a model with freely varying e parameter
summary(model1v)


rm(model1ca, model1h, model1v)
}
\keyword{datasets}
