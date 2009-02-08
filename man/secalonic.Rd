\name{secalonic}

\alias{secalonic}

\docType{data}

\title{Root length}

\description{
  Data are from an experiment assessing the inhibitory effect of secalonic acids on plant growth.
}

\usage{data(secalonic)}

\format{
  A data frame with 7 observations on the following 2 variables:
  \describe{
    \item{\code{dose}}{a numeric vector containing dose values}
    \item{\code{rootl}}{a numeric vector containing average root lengths}
  }
}

\details{
  Data are averages of three root length measurements.
}

\source{
  Gong, X. and Zeng, R. and Luo, S. and Yong, C. and Zheng, Q. (2004) Two new 
  secalonic acids from \emph{Aspergillus Japonicus} and their allelopathic effects on higher plants,
  \emph{Proceedings of International Symposium on Allelopathy Research and Application, 27-29 April, 
  Shanshui, Guangdong, China (Editors: R. Zeng and S. Luo)}, 209--217.
}

%\references{}

\examples{

## Fitting a 4-parameter log-logistic model
secalonic.m1 <- multdrc(rootl~dose, data=secalonic)

## Fitting a 3-parameter log-logistic model
##  lower limit fixed at 0
secalonic.m2 <- multdrc(rootl~dose, data=secalonic, fct=LL.3())


## Fitting logistic and log-logistic models
secalonic.m1 <- drm(rootl~dose, data=secalonic, fct=LL.4())
secalonic.m2 <- drm(rootl~dose, data=secalonic, fct=L.4())

plot(secalonic.m1, conLevel=0.001, broken=TRUE, ylim=c(0,7), xlab="Dose", ylab="Root length", cex=1.2, lwd=2)
plot(secalonic.m2, conLevel=0.001, broken=TRUE, ylim=c(0,7), add=TRUE, lty=2, lwd=2)

abline(h=coef(secalonic.m2)[3], lty=3, lwd=2)

}
\keyword{datasets}
