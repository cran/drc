\name{spinach}

\alias{spinach}

\docType{data}

\title{Inhibition of photosynthesis}

\description{
  Experiment investigating the inhibition of photosynthesis in response to two synthetic 
  photosystem II inhibitors, the herbicides diuron and bentazone.
  The effect of oxygen consumption of thylakoid membranes (chloroplasts) from spinach
  was measured after incubation with the synthetic inhibitors, resulting in 5 dose-response curves.
}

\usage{data(spinach)}

\format{
  A data frame with 105 observations on the following 3 variables:
  \describe{
    \item{CURVE}{a numeric vector.}
    \item{HERBICIDE}{a character vector.}
    \item{DOSE}{a numeric vector.}
    \item{SLOPE}{a numeric vector.}
  }
}

\details{
  CURVE denotes a total of five independent assays with two herbicides.
  HERBICIDE denotes the names of the herbicides.
  DOSE is the concentration of the herbicides muMol and SLOPE is 
  the oxygen consumption of thylakoid membranes.  
}

\source{
  Streibig, J. C. (1998) Joint action of natural and synthetic photosystem II inhibitors, \emph{Pesticide Science}, \bold{55}, 137--146.
}

\examples{

## Displaying the first rows in the dataset
head(spinach) [1:5,]  # displaying first 5 rows in the data set

}
\keyword{datasets}
