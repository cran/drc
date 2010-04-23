\name{BIC-methods}

\docType{methods}

\alias{BIC-methods}

\alias{BIC,drc-method}

\title{
  Bayesian Information Criterion for dose-response models
}

\description{
  This method calculates Schwarz's Bayesian criterion for dose-response model fits.
}

\section{Methods}{
\describe{

\item{object = "drc"}{Bayesian Information Criterion for dose-response models}
}}

\author{
  Tobias Verbeke and Christian Ritz
}

%\note{}

\seealso{
  The help page for \code{\link[nlme]{BIC}} provides more details.
}

\examples{

CONC <- c(6,6.52,7,7.52,8,8.52,9,10)
RESP <- c(15.01,21.33,19.58,48.35,70.86,75.3,80.16,81.7)

nlsEx <- nls(RESP ~ SSfpl(CONC, A, B, xmid, scale))
BIC(nlsEx)
# [1] 46.9176

drmEx <- drm(RESP ~ CONC, fct = L.4())
BIC(drmEx)
# [1] 46.9176


}
\keyword{methods}
\keyword{models}
\keyword{nonlinear} 
