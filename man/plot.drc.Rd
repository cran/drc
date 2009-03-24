\name{plot.drc}

\alias{plot.drc}
\alias{plot.mrdrc}

\title{Plotting fitted curves for a 'drc' or 'mrdrc' object}

\description{
  \code{plot} displays fitted curves and observations in the same plot window, 
  distinguishing between curves by different plot symbols and line types.
}

\usage{
  
  \method{plot}{drc}(x, ..., level = NULL, broken = FALSE, col = FALSE, 
  conLevel, conName, gridsize = 100, legend, legendText, legendPos, legendCex = 1,
  type = c("average", "all", "bars", "none", "obs"), obs, lty, log = "x", xsty = c("base10", "none"), 
  cex, pch, xlab, ylab, xlim, ylim, cex.axis = 1,
  bcontrol = NULL, xt = NULL, xtlab = NULL, yt = NULL, ytlab = NULL,
  add = FALSE, axes = TRUE)

  \method{plot}{mrdrc}(x, ..., pava = FALSE)  
}

\arguments{
  \item{x}{an object of class 'drc'.}
  \item{...}{additional graphical arguments. For instance, use \code{lwd=2} or \code{lwd=3} to increase the width 
  of plot symbols.}
  \item{level}{vector of character strings. To plot only the curves specified by their names.}
  \item{broken}{logical. If TRUE the x axis is broken provided this axis is logarithmic 
  (using functionality in the CRAN package 'plotrix').}
  \item{col}{either logical or a vector of colours. If TRUE default colours are used. If FALSE (default) no colours are used.}  
  \item{conLevel}{numeric specifying the dose level below which the dose is zero (the amount of stretching on 
  the dose axis above zero in order to create the visual illusion of a logarithmic scale including 0). 
  The default is the base-10 value corresponding to the rounded value of the minimum of the log10 values of 
  all positive dose values.}
  \item{conName}{character string. Name on x axis for dose zero. Default is '"0"'.}
  \item{gridsize}{numeric. Number of points in the grid used for plotting the fitted curves.}
  \item{legend}{logical. If TRUE a legend is displayed.}
  \item{legendText}{a character string or vector of character strings specifying the legend text 
  (the position of the upper right corner of the legend box).}
  \item{legendPos}{numeric vector of length 2 giving the position of the legend.}
  \item{legendCex}{numeric specifying the legend text size.}
  \item{type}{a character string specifying how to plot the data. There are currently
  5 options:  "average" (averages and fitted curve(s); default), 
  "none" (only the fitted curve(s)), "obs" (only the data points), 
  "all" (all data points and fitted curve(s)) and 
  "bars" (averages and fitted curve(s) with model-based error ).}
  \item{obs}{Outdated argument. Use \code{type}.}  
  \item{lty}{a numeric vector specifying the line types.}
  \item{log}{a character string which contains '"x"' if the x axis is to be logarithmic, '"y"' if the y axis is to be logarithmic and '"xy"' or
  '"yx"' if both axes are to be logarithmic. The default is "x". The empty string "" yields the original axes.}
  \item{xsty}{a character string specifying the dose axis style for arrangement of tick marks. By default ("base10") 
  only base 10 tick marks are shown, unless too few tick marks are available in which case "none" applies 
  (if the range of doses applied is narrow). The choice "none" seeks to provide a sensible set of tick marks.}
  \item{cex}{numeric or numeric vector specifying the size of plotting symbols and text 
  (see \code{\link{par}} for details).}
  \item{pch}{a vector of plotting characters or symbols (see \code{\link{points}}).}
  \item{xlab}{an optional label for the x axis.}
  \item{ylab}{an optional label for the y axis.}
  \item{xlim}{a numeric vector of length two, containing the lower and upper limit for the x axis.}
  \item{ylim}{a numeric vector of length two, containing the lower and upper limit for the y axis.}
  \item{cex.axis}{numeric value specifying the magnification to be used for axis annotation 
  relative to the current setting of cex.}  
  \item{bcontrol}{a list with components \code{factor}, \code{style} and \code{width}. 
  Controlling the appearance of the break (in case \code{broken} is \code{TRUE}). 
  The component \code{factor} is the distance from the control to the break as a 
  multiple of the value of \code{conLevel} (default is 2).
  The component \code{style} can take the values: \code{gap}, \code{slash} or \code{zigzag}.
  The component \code{width} is the width of the break symbol (default is 0.02).}
  \item{xt}{a numeric vector containing the positions of the tick marks on the x axis.}
  \item{xtlab}{a vector containing the tick marks on the x axis.}
  \item{yt}{a numeric vector, containing the positions of the tick marks on the y axis.}
  \item{ytlab}{a vector containing the tick marks on the y axis.} 
  \item{add}{logical. If TRUE then add to already existing plot.} 
  \item{axes}{logical indicating whether both axes should be drawn on the plot.}
  \item{pava}{logical. If TRUE the fit is monotoniosed using pool adjacent violators algorithm.}  
}

\details{
   \code{plot.mrdrc} takes the same arguments as \code{plot.drc}.

  Suitable labels are automatically provided.
  
  The use of \code{xlim} allows changing the range of the x axis, extrapolating the fitted dose-response curves. 
  Note that changing the range on the x axis may also entail a change of the range on the y axis. Sometimes
  it may be useful to extend the upper limit on the y axis (using \code{ylim}) in order to fit a legend into 
  the plot.
  
  See \code{\link{colors}} for the available colours.
  
  The arguments \code{broken} and \code{bcontrol} rely on the function \code{link{axis.break}} with arguments 
  \code{style} and \code{brw} in the package \code{plotrix}.
}

\value{
  An invisible data frame with the values used for plotting the fitted curves. The first column contains the dose values,
  and the following columns (one for each curve) contain the fitted response values.
}

\author{
  Christian Ritz and Jens C. Streibig
  Contributions from Xiaoyan Wang  
}

%\note{ ~~further notes~~ }

%\seealso{\code{\link{plotraw}} plots the observations only.}

\examples{

## Fitting models to be plotted below
ryegrass.m1 <- drm(rootl~conc, data = ryegrass, fct = LL.4())
ryegrass.m2 <- drm(rootl~conc, data = ryegrass, fct = LL.3())  # lower limit fixed at 0

## Plotting observations and fitted curve for the first model
plot(ryegrass.m1, broken = TRUE) 

## Adding fitted curve for the second model (not much difference)
plot(ryegrass.m2, broken = TRUE, add = TRUE, type = "none", col = 2, lty = 2)

## Finetuning the axis break
plot(ryegrass.m1, broken = TRUE, bcontrol = list(style = "gap"))
plot(ryegrass.m1, broken = TRUE, bcontrol = list(style = "slash"))
plot(ryegrass.m1, broken = TRUE, bcontrol = list(style = "zigzag"))

## Plot without axes
plot(ryegrass.m1, axes = FALSE)

## Fitting model to be plotted below
spinach.m1 <- drm(SLOPE~DOSE, CURVE, data = spinach, fct = LL.4())

## Plot with no colours
plot(spinach.m1, main = "Different line types (default)")

## Plot with default colours
plot(spinach.m1, col = TRUE, main = "Default colours")

## Plot with specified colours
plot(spinach.m1, col = c(2,6,3,23,56), main = "User-specified colours")

## Plot of curves 1 and 2 only
plot(spinach.m1, level = c(1,2), main = "User-specified curves")

## Plot with symbol of different sizes
plot(spinach.m1, cex = c(1,2,3,4,5), main = "User-specified symbil sizes")

## Fitting another model to be plotted below
lettuce.m1 <- drm(weight~conc, data = lettuce, fct = LL.4())

## Using the argument 'conLevel'. Compare the plots!
par(mfrow = c(2, 2))
plot(lettuce.m1, main = "conLevel=1e-2 (default)")  # using the default
plot(lettuce.m1, conLevel = 1e-4, main = "conLevel=1e-4")
plot(lettuce.m1, conLevel = 1e-6, main = "conLevel=1e-6")
plot(lettuce.m1, conLevel = 1e-8, main = "conLevel=1e-8")
par(mfrow = c(1,1))

## User-specified position of legend
S.alba.m1 <- drm(DryMatter~Dose, Herbicide, data = S.alba, fct = LL.4())

plot(S.alba.m1)
plot(S.alba.m1, legendPos = c(0.3, 4.8))

}
\keyword{aplot}
