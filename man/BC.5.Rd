\name{BC.5}

\alias{BC.5}
\alias{bcl4}

\alias{BC.4}
\alias{bcl3}

\title{The Brain-Cousens hormesis models}

\description{
  'BC.4' and 'BC.5' provide the Brain-Cousens modified log-logistic models for describing u-shaped hormesis.
}

\usage{
  BC.5(names = c("b", "c", "d", "e", "f"))

  BC.4(names = c("b", "d", "e", "f"))
}

\arguments{
  \item{names}{a vector of character strings giving the names of the parameters. The default is reasonable.}
}

\details{
  The Brain-Cousens model is given by the expression
  \deqn{ f(x) = c + \frac{d-c+fx}{1+\exp(b(\log(x)-\log(e)))}}
  which is a five-parameter model.
  
  It is a modification of the log-logistic curve to take u-shaped hormesis into account.
  
  Fixing the lower limit at 0 yields the four-parameter model 
  \deqn{ f(x) = 0 + \frac{d-0+fx}{1+\exp(b(\log(x)-\log(e)))}}
}

\value{
  See \code{\link{braincousens}}.
}

\references{
  Brain, P. and Cousens, R. (1989) An equation to describe dose responses 
  where there is stimulation of growth at low doses, 
  \emph{Weed Research}, \bold{29}, 93--96.

  van Ewijk, P. H. and Hoekstra, J. A. (1993) 
  Calculation of the EC50 and its Confidence Interval When Subtoxic Stimulus Is Present,
  \emph{ECOTOXICOLOGY AND ENVIRONMENTAL SAFETY}, \bold{25}, 25--32.  
}

\author{Christian Ritz}

\note{This function is for use with the function \code{\link{drm}} or \code{\link{multdrc}}.}

%\seealso{The related function \code{\link{bcl3}}.}

\examples{

model1 <- multdrc(hormesis[,c(2,1)], fct=BC.5())
anova(model1)
plot(model1)

model2 <- multdrc(hormesis[,c(2,1)], fct=BC.4())
summary(model2)
ED(model2, c(50))  # compare the parameter estimate and 
                   # its estimated standard error 
                   # to the values in the paper by 
                   # van Ewijk and Hoekstra

rm(model1, model2)

}
\keyword{models}
\keyword{nonlinear}
