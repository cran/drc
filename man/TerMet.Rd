\name{TerMet}

\alias{TerMet}

\docType{data}

\title{The effect of terbuthylazin on growth rate}

\description{
  Test on the effect of terbuthylazin on  \emph{Lemna minor}, performed on an aseptic
  culture according to the OECD-guidelines.  
}

\usage{data(TerMet)}

\format{
  A data frame with 30 observations on the following 2 variables.
  \describe{
    \item{dose}{a numeric vector}
    \item{rgr}{a numeric vector}
  }
}

\details{
  Dose is \deqn{\mu l^{-1}} and rgr is the relative growth rate of \emph{Lemna}.  
}

\source{
  Cedergreen N. (2004) Unpublished bioassay data.
}

\examples{

## displaying first 6 rows of the data set
head(TerMet)

## Fitting log-logistic model
m1 <- drm(rgr~dose, data = TerMet, fct = LL.4())
summary(m1)

## Fitting log-logistic model
##  with Box-Cox transformation
m1 <- drm(rgr~dose, data = TerMet, fct = LL.4(), adjust = "bc1")
summary(m1)

}
\keyword{datasets}
