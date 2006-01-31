\name{plot.drc}

\alias{plot.drc}

\title{Plotting fitted curves for a drc object}

\description{
  'plot' displays fitted curves and observations in the same plot window, distinguishing between curves by different plot symbols and line types.
}

\usage{
  plot.drc(x, ..., level=NULL, breakCurve = FALSE, col=FALSE, 
  conLevel, conName, grid=100, legend=TRUE, legendText, 
  type = "average", obs, lty, log="x", pch, xlab, ylab, xlim, ylim)
}

\arguments{
  \item{x}{an object of class 'drc'. For instance, 'lwd=2' or 'lwd=3' increase the width of plot symbols.}
  \item{...}{additional arguments.}
  \item{level}{vector of character strings. To plot only the curves specified by their names.}
  \item{breakCurve}{logical. If TRUE the curves break before reaching control. Experimental.}
  \item{col}{either logical or a vector of colours. If TRUE default colours are used. If FALSE (default) no colours are used.}  
  \item{conLevel}{numeric. Dose level below which the dose is zero (the amount of stretching on the x-axis above zero). Default is 1e-2.}
  \item{conName}{character string. Name on x axis for dose zero. Default is '"0"'.}
  \item{grid}{numeric. Number of points in the grid used for plotting the fitted curves.}
  \item{legend}{logical. If TRUE a legend is displayed.}
  \item{legendText}{a character string or vector of character strings specifying the legend text.}
  \item{type}{a character string specifying how the originals observations should be plotted. There are 4 options: "average" (default), "none", "points"
  and "add" (for adding a fit to an existing plot).}
  \item{obs}{Outdated argument. Use \code{type}.}  
  \item{lty}{a numeric vector specifying the line types.}
  \item{log}{a character string which contains '"x"' if the x axis is to be logarithmic, '"y"' if the y axis is to be logarithmic and '"xy"' or
  '"yx"' if both axes are to be logarithmic. The default is "x". The empty string "" yields the original axes.}
  \item{pch}{a vector of plotting characters or symbols (see \code{\link{points}}).}
  \item{xlab}{an optional label for the x axis.}
  \item{ylab}{an optional label for the y axis.}
  \item{xlim}{a numeric vector of length two, containing the lower and upper limit for the x axis.}
  \item{ylim}{a numeric vector of length two, containing the lower and upper limit for the y axis.}  
}

\details{
  Suitable labels are automatically provided.
  
  The use of \code{xlim} allows changing the x-axis, extrapolating the fitted dose-response curves.
  
  See \code{\link{colors}} for the available colours.
}

\value{
  An invisible data frame with the values used for plotting the fitted curves. The first column contains the dose values,
  and the following columns (one for each curve) contain the fitted response values.
}

\author{Christian Ritz and Jens C. Streibig}

%\note{ ~~further notes~~ }

\seealso{\code{\link{plotraw}} plots the observations only.}

\examples{

## Fitting models to be plotted below
model1 <- multdrc(MEANLR~MM, data=FA)
model2 <- multdrc(MEANLR~MM, data=FA, fct=l3())  # lower limit fixed at 0

## Plotting observations and fitted curve for 'model1'
plot(model1)

## Adding fitted curve for 'model2'
plot(model2, type="add")

## Fitting model to be plotted below
model3 <- multdrc(SLOPE~DOSE, CURVE, data=PestSci)

## Plot with no colours
plot(model3)

## Plot with default colours
plot(model3, colour=TRUE)

## Plot with specified colours
plot(model3, colour=c(2,6,15,23,56))

## Plot of curves 1 and 2 only
plot(model3, level=c(1,2))

## Fitting another model to be plotted below
model4 <- multdrc(weight~conc, data=hormesis)

## Using the argument 'conLevel'. Compare the plots!
plot(model4)  # using default setting 'conLevel=1e-2'
plot(model4, conLevel=1e-4)
plot(model4, conLevel=1e-8)

## Removing models from work space
rm(model1, model2, model3, model4)

}
\keyword{aplot}
