\name{plotraw}

\alias{plotraw}

\title{Plotting raw dose response data.}

\description{
  'plotraw' plots raw dose response data, either all curves in a single coordinate system or in separate coordinate systems.
}

\usage{
  plotraw(formula, curveid, data = NULL, na.action = na.fail, trtid, 
  colour = FALSE, conLevel = 0.01, conName = "0", trellis = FALSE, 
  log = "x", xlab, ylab, ...)
}

\arguments{
  \item{formula}{a symbolic description of the model to be fit. Either of the form 'response \eqn{~} dose' or as a data frame with response value
                 in first column and dose in second column.}
  \item{curveid}{a vector containing the grouping of data into separate curves.}
  \item{data}{an optional data frame where dose, response and curveid are found.}
  \item{na.action}{a function which indicates what should happen when the data contain 'NA's.  The default is 'na.fail'. To omit 'NA's use
  'na.omit'.}
  \item{trtid}{a character vector. The names of the treatment underlying the curves.}
  \item{colour}{either logical or a vector of colours. If TRUE default colours are used. If FALSE (default) no colours are used.}  
  \item{conLevel}{numeric. Dose level below which the dose is zero (the amount of stretching on the x-axis above zero). Default is 1e-2.}
  \item{conName}{character. Name on x axis for dose zero. Default is '"0"'.}
  \item{trellis}{logical. If TRUE a Trellis plot is produced. If FALSE all curves are displayed in the same coordinate system.}
  \item{log}{a character string which contains '"x"' if the x axis is to be logarithmic, '"y"' if the y axis is to be logarithmic and '"xy"' or
  '"yx"' if both axes are to be logarithmic. The default is "x". The empty string "" yields the original axes.}
  \item{xlab}{an optional label for the x axis.}
  \item{ylab}{an optional label for the y axis.}
  \item{...}{additional arguments.}
}
%\details{}

%\value{}

%\references{}

\author{Christian Ritz}

%\note{ ~~further notes~~ }

\seealso{See Also as \code{\link{plot.drc}}.}

%\examples{}

\keyword{aplot}
