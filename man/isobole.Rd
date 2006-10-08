\name{isobole}
\alias{isobole}

\title{Creating isobolograms}

\description{
  'isobole' displays isobole based on EC/ED50 estimates from a log-logistic model.
  Additionally isoboles determined by the concentration addition model or Hewlett's
  model can be added to the plot.
}

\usage{
isobole(object1, object2, ename = "e", xaxis = "0", exchange, xlab = "A", ylab = "B", xlim, ylim, ...)
}

%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{object1}{object of class 'drc' where EC/ED50 parameters vary freely.}
  \item{object2}{object of class 'drc' where EC/ED50 parameters vary according to Hewlett's model.}
  \item{ename}{character string. The name of the EC/ED50 variable.}
  \item{xaxis}{character string. Is the mixture "0:100" or "100:0" on the x axis?}
  \item{exchange}{numeric. The exchange rate between the two substances.}
  \item{xlab}{an optional label for the x axis.}
  \item{ylab}{an optional label for the y axis.}
  \item{xlim}{a numeric vector of length two, containing the lower and upper limit for the x axis.}
  \item{ylim}{a numeric vector of length two, containing the lower and upper limit for the y axis.}
  \item{\dots}{Additional graphical parameters.}
}

\details{
  Only to be used with the function \code{\link{mixture}}.
}

%\value{}

%\references{ ~put references to the literature/web site here ~ }

\author{Christian Ritz}

%\note{}

%\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }

\examples{

## Fitting the model with freely varying ED50 values
model1 <- multdrc(rgr1~dose1, pct1, data=acidiq, 
collapse = list(~factor(pct1),~1,~1,~factor(pct1) - 1), 
boxcox = TRUE)

## Fitting the concentration addition model
model1ca <- mixture(rgr1~dose1, pct1, data=acidiq, 
collapse = ~factor(pct1), 
model = "CA", boxcox = TRUE)

## Fitting the Hewlett model
model1h <- mixture(rgr1~dose1, pct1, data=acidiq, 
collapse = ~factor(pct1), 
model = "Hewlett", boxcox = TRUE)

## Isobologram with freely varying EC/ED50 values
isobole(model1, exchange=0.1, xlim=c(0,450), ylim=c(0,45))

## Isobologram based on concentration addition
isobole(model1, model1ca, exchange=0.1, xlim=c(0,450), ylim=c(0,45), 
main="Concentration Addition")

## Isobologram based on Hewlett's model
isobole(model1, model1h, exchange=0.1, xlim=c(0,450), ylim=c(0,45), 
main="Hewlett")

rm(model1, model1ca, model1h)

}
\keyword{models}
\keyword{nonlinear}
