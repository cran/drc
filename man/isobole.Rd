\name{isobole}
\alias{isobole}

\title{Creating isobolograms}

\description{
  'isobole' displays isobole based on EC/ED50 estimates from a log-logistic model.
  Additionally isoboles determined by the concentration addition model, Hewlett's
  model and Voelund's model can be added to the plot.
}

\usage{
isobole(object1, object2, ename = "e", xaxis = "0", exchange = 1, cifactor = 2, xlab = "0", ylab = "100", xlim, ylim, ...)
}

\arguments{
  \item{object1}{object of class 'drc' where EC/ED50 parameters vary freely.}
  \item{object2}{object of class 'drc' where EC/ED50 parameters vary according to Hewlett's model.}
  \item{ename}{character string. The name of the EC/ED50 variable.}
  \item{xaxis}{character string. Is the mixture "0:100" or "100:0" on the x axis?}
  \item{exchange}{numeric. The exchange rate between the two substances.}
  \item{cifactor}{numeric. The factor to be used in the confidence intervals. 
  Default is 2, but 1 has been used in publications.}  
  \item{xlab}{an optional label for the x axis.}
  \item{ylab}{an optional label for the y axis.}
  \item{xlim}{a numeric vector of length two, containing the lower and upper limit for the x axis.}
  \item{ylim}{a numeric vector of length two, containing the lower and upper limit for the y axis.}
  \item{\dots}{Additional graphical parameters.}
}

\details{
  The model fits to be supplied as first and optionally second argument are obtained
  using \code{\link{mixture}} and \code{\link{drm}}.
}

\value{
  No value is returned. Only used for the side effect: the isobologram shown.
}

\references{
  Soerensen, H. and Cedergreen, N. and Skovgaard, I. M. and Streibig, J. C. (2007) 
  An isobole-based statistical model and test for synergism/antagonism in binary mixture toxicity experiments, 
  \emph{Environmental and Ecological Statistics}, \bold{14}, 383--397.
}

\author{Christian Ritz}

%\note{}

\seealso{
  The examples in \code{\link{acidiq}} (the Hewlett model), \code{\link{glymet}} (concentration addition) 
  and \code{\link{mecter}} (the Voelund model).
}

\examples{
}
\keyword{models}
\keyword{nonlinear}
