\name{christensenetal2003}

\alias{christensenetal2003}

\docType{data}

\title{Potency of two herbicides}

\description{
  Data are from an experiment, comparing the potency of the two herbicides glyphosate and bentazone in
  \emph{Sinapis alba}.
}

\usage{data(christensenetal2003)}

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
  The dataset is analysed in Christensen \emph{et al} (2003).
}

\references{
  Christensen, M. G. and Teicher, H. B. and Streibig, J. C. (2003) Linking fluorescence 
  induction curve and biomass in herbicide screening, \emph{Pest Management Science},
  \bold{59},  1303--1310.
}

\examples{

## Fitting a log-logistic model with
##  common lower and upper limits
model1<-multdrc(DryMatter~Dose, Herbicide, data=christensenetal2003,
collapse=data.frame(Herbicide,1,1,Herbicide), boxcox=TRUE) 

## Plotting fitted curve and data
plot(model1)


rm(model1)
}
\keyword{datasets}
