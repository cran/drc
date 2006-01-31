\name{onecurve}

\alias{onecurve}

\docType{data}

\title{Root length}

\description{
  Data are from an experiment assessing the inhibitory effect effect of secalonic acids on plant growth.
}

\usage{data(onecurve)}

\format{
  A data frame with 7 observations on the following 2 variables.
  \describe{
    \item{\code{Dose}}{a numeric vector containing dose values}
    \item{\code{Root}}{a numeric vector containing root lengths}
  }
}

\details{
  The root lengths are each averages of three measured root lengths.
}

\source{
  The dataset is analysed in Gong \emph{et al} (2004).
}

\references{
  Gong, X. and Zeng, R. and Luo, S. and Yong, C. and Zheng, Q. (2004) Two new 
  secalonic acids from \emph{Aspergillus Japonicus} and their allelopathic effects on higher plants,
  \emph{Proceedings of International Symposium on Allelopathy Research and Application, 27-29 April, 
  Shanshui, Guangdong, China (Editors: R. Zeng and S. Luo)}, 209--217.
}

\examples{

## Fitting a 4-parameter log-logistic model
model.onecurve <- multdrc(Root~Dose, data=onecurve)

## Fitting a 3-parameter log-logistic model
##  lower limit fixed at 0
model.onecurve.noC <- multdrc(Root~Dose, data=onecurve, fct=l3())


rm(model.onecurve, model.onecurve.noC)
}
\keyword{datasets}
