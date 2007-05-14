\name{cedergreen}

\alias{cedergreen}
\alias{CRS.6}
\alias{ucedergreen}

\title{The Cedergreen-Ritz-Streibig model}

\description{
  'cedergreen' provides a very general way of specifying then Cedergreen-Ritz-Streibig 
  modified log-logistic model for describing hormesis, under various constraints on the parameters.

  \code{\link{CRS.6}} is the extension of \code{link{cedergreen}} with freely varying alpha parameter.

  For u-shaped hormesis data 'ucedergreen' provides a very general way of specifying the 
  Cedergreen-Ritz-Streibig modified log-logistic model, under various constraints on the parameters.
}

\usage{
  cedergreen(lowerc = c(-Inf, -Inf, -Inf, -Inf, -Inf), 
  upperc = c(Inf, Inf, Inf, Inf, Inf), fixed = c(NA, NA, NA, NA, NA), 
  names = c("b", "c", "d", "e", "f"), alpha, scaleDose = TRUE)
  
  CRS.6(lowerc=c(-Inf, -Inf, -Inf, -Inf, -Inf, -Inf), 
  upperc=c(Inf, Inf, Inf, Inf, Inf, Inf), 
  fixed=c(NA, NA, NA, NA, NA, NA), 
  names=c("b","c","d","e","f","g"))
  
  ucedergreen(lowerc = c(-Inf, -Inf, -Inf, -Inf, -Inf), 
  upperc = c(Inf, Inf, Inf, Inf, Inf), fixed = c(NA, NA, NA, NA, NA), 
  names = c("b", "c", "d", "e", "f"), alpha, scaleDose = TRUE)
}

\arguments{
  \item{lowerc}{numeric vector. The lower bound on parameters. Default is minus infinity.}
  \item{upperc}{numeric vector. The upper bound on parameters. Default is plus infinity.}
  \item{fixed}{numeric vector. Specifies which parameters are fixed and at what value they are fixed. 
  NAs for parameter that are not fixed.}
  \item{names}{a vector of character strings giving the names of the parameters (should not contain ":"). 
  The default is reasonable (see under 'Usage').
               The order of the parameters is: b, c, d, e, f (see under 'Details').}
  \item{alpha}{numeric. The degree of hormesis. Needs to be specified!}
  \item{scaleDose}{logical. If TRUE dose values are scaled around 1 during estimation; 
                   this is required for datasets where all dose values are small. }  
}

\details{
The model is given by the expression
  
  \deqn{ f(x) = c + \frac{d-c+f exp(-1/(x^{\alpha}))}{1+exp(b(log(x)-log(e)))}}
  
which is a five-parameter model (alpha is fixed or freely varying). Not all features (eg EC/ED calculation) 
are available for the model with freely varying alpha.
  
It is a modification of the four-parameter logistic curve to take hormesis into account.
  
The u-shaped model is given by the expression
  
  \deqn{ f(x) = cd - \frac{d-c+f \exp(-1/x^{\alpha})}{1+\exp(b(\log(x)-\log(e)))}}
  
}

\value{
  The value returned is a list containing the non-linear function, the self starter function
  and the parameter names.
}

\references{
  Cedergreen, N. and Ritz, C. and Streibig, J. C. (2005) Improved empirical models describing hormesis, 
  \emph{Environmental Toxicology and Chemistry} \bold{24}, 3166--3172.
}

\author{Christian Ritz}

\note{The functions are for use with the functions \code{\link{drm}} or \code{\link{multdrc}}.}

\seealso{
  Special cases are \code{\link{CRS.4a}}, \code{\link{CRS.4a}}, 
  \code{\link{UCRS.5a}} and \code{\link{UCRS.5a}} where a,b and c coresspond to 
  the pre-specified alpha values 1, 0.5 and 0.25, respectively.
}

\examples{

## Modified logistic model with the constraint f>0
model1 <- multdrc(hormesis[,c(2,1)], 
fct=cedergreen(fixed=c(NA, NA, NA, NA, NA), 
lowerc=c(-Inf, -Inf, -Inf, -Inf, 0), alpha=1), 
control=mdControl(constr=TRUE))

summary(model1)
ED(model1, c(10, 50, 90))


## Estimating alpha
m2<-drm(weight~conc,data=hormesis,fct=CRS.6())
summary(m2)
plot(m2)  # oops: not increasing until hormesis peak

}

\keyword{models}
\keyword{nonlinear}

\concept{hormesis hormetic effect initial stimulation u-shaped}
