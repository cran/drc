\name{PestSci}

\alias{PestSci}

\docType{data}

\title{Inhibition of photosynthesis}

\description{
  Inhibition of photosynthesis in response to two synthetic photosystem II inhibitors, the herbicides diuron and bentazone.
  The effect of oxygen evolution of thylakoid membranes (Chloroplasts) from Spinach
  was measured after incubation with the synthetic inhibitors, resulting in 5 dose-response curves.
}

\usage{data(PestSci)}

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
  the Oxygen evolution from thylakoid membranes.  
}

\source{
  Streibig, J. C. (1998) Joint action of natural and synthetic photosystem II inhibitors, \emph{Pesticide Science}, \bold{55}, 137--146.
}

\examples{

PestSci[1:5,]  # displaying first 5 rows in the data set
str(PestSci)   # compactly displaying the structure of the PestSci object

## Analysis accounting for variance heterogeneity 
##  using Box-Cox-transformation
model1 <- multdrc(SLOPE~DOSE, CURVE, data=PestSci, boxcox=TRUE)
summary(model1)  
## note: some estimates have indeterminate estimated standard error

model2 <- multdrc(SLOPE~DOSE, CURVE, data=PestSci, boxcox=TRUE,
control=mdControl(bcAdd=0.05))
summary(model2)  
## note: now all estimates have estimated standard errors

## Analysis accounting for variance heterogeneity using variance modelling
model3 <- multdrc(SLOPE~DOSE, CURVE, data=PestSci, varPower=TRUE)
summary(model3)

## The same model in specified in 4 ways (first 3 with same parameterisation)
model4 <- multdrc(SLOPE~DOSE, CURVE, data=PestSci)

model5 <- multdrc(SLOPE~DOSE, CURVE, data=PestSci, collapse=data.frame(CURVE,CURVE,CURVE,CURVE))

model6 <- multdrc(SLOPE~DOSE, CURVE, data=PestSci, 
collapse=list(~factor(CURVE)-1,~factor(CURVE)-1,~factor(CURVE)-1,~factor(CURVE)-1))

model7 <- multdrc(SLOPE~DOSE, CURVE, data=PestSci, 
collapse=list(~factor(CURVE),~factor(CURVE),~factor(CURVE),~factor(CURVE)))

## Example using 'hetvar' argument 
##  for heterogeneous variances
model8 <- multdrc(SLOPE~DOSE, CURVE, data = PestSci, 
hetvar = HERBICIDE, boxcox = TRUE, bcAdd = 0.01)
summary(model8)

rm(model1, model2, model3, model4, model5, model6, model7, model8)

}
\keyword{datasets}
