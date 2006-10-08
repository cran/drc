\name{leaflength}

\alias{leaflength}

\docType{data}

\title{Leaf length of barley}

\description{
  In an experiment barley was grown in a hydroponic solution with a herbicide.
}

\usage{data(leaflength)}

\format{
  A data frame with 42 observations on the following 2 variables.
  \describe{
    \item{\code{Dose}}{a numeric vector}
    \item{\code{DW}}{a numeric vector}
  }
}

\details{
  The dataset exhibits a large hormetical effect.
}

\source{
  Nina Cedergreen, Royal Veterinary and Agricultural University, Denmark.
}

%\references{}

\examples{

## Fitting a hormesis model
fit.ml3c <- multdrc(DW~Dose, data=leaflength, fct=ml3c(), boxcox=TRUE)

## Plottinf fitted curve and original data
plot(fit.ml3c, conLevel=0.001, type="points", legend=FALSE, 
ylab="Produced leaf length(cm)", xlab="Metsulfuron-methyl(mg/l)",
main="Hormesis: leaf length of barley")

rm(fit.ml3c)

}
\keyword{datasets}
