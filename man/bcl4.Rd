\name{bcl4}
\alias{bcl4}

\title{The four-parameter logistic function with Brain-Cousens' modification}

\description{
  'bcl4' provides the five-parameter function, self starter function and names of the parameters.
}

\usage{
bcl4(names = c("b", "c", "d", "e", "f"), useDer = FALSE)
}

\arguments{
  \item{names}{a vector of character strings giving the names of the parameters. The default is reasonable.}
  \item{useDer}{logical. If TRUE derivatives are supplied, otherwise they are not supplied. Not implemented!}
}

\details{
  The Brain-Cousens model is given by the expression
  \deqn{ f(x) = c + \frac{d-c+fx}{1+\exp(b(\log(x)-\log(e)))}}
  which is a five-parameter model.
  
  It is a modification of the four-parameter logistic curve to take hormesis into account.
}

\value{
  See \code{\link{braincousens}}.
}

\references{Brain, P. and Cousens, R. (1989) An equation to describe dose responses where there is stimulation of growth at low doses, 
            \emph{Weed Research}, \bold{29}, 93--96.}

\author{Christian Ritz}

\note{This function is for use with the function \code{\link{multdrc}}.}

\seealso{The related function \code{\link{bcl3}}.}

\examples{

model1 <- multdrc(hormesis[,c(2,1)], fct=bcl4())
anova(model1)
plot(model1)

rm(model1)

}
\keyword{models}
\keyword{nonlinear}
