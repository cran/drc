\name{MM}
\alias{MM}

\title{Michaelis-Menten model}

\description{
  'MM' specifies the Michaelis-Menten model or the shifted Michaelis-Menten model
  for modelling enzyme kinetics, weed densities etc.
}

\usage{
MM(fixed = c(NA, NA, 0), names = c("Vm", "K", "y0"))
}

\arguments{
  \item{fixed}{numeric vector. Specifies which parameters are fixed and at what value they are fixed. 
  NAs for parameter that are not fixed.}
  \item{names}{a vector of character strings giving the names of the parameters (should not contain ":"). 
  The order of the parameters is: Vm, K, y0 (see under 'Details').}}

\details{
  The model is given by the expression
  
  \deqn{ f(x) = y0 + \frac{Vm x}{K + x}}
  
  It is a three-parameter model.
}

\value{
  An object containing with all components needed in \code{\link{multdrc}}.
}

%\references{}

\author{Christian Ritz}

\note{
  At the moment the implementation cannot deal with infinite concentrations.
}

%\seealso{}

\examples{

ctl <- data.frame(list(conc=c(0,0,0), rate=c(55,60,65), state=c('ctl','ctl','ctl')))
pu <- rbind(Puromycin, ctl)

model1<-multdrc(rate~conc, state, 
collapse=data.frame(state,state,1), 
data=pu, fct=MM(fixed=c(NA,NA,NA)))

summary(model1)
plot(model1) 

rm(ctl, pu, model1)

}
\keyword{models}
\keyword{nonlinear}
