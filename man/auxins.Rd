\name{auxins}

\alias{auxins}

\docType{data}

\title{ Effect of technical grade and commerecially formulated auxin herbicides}

\description{
 MCPA, 2,4-D, mecorprop and dichorlprop were applied either as technical grades
 materials  (h = 1, 2, 3, 4) or as commercial formulations (herb = 5, 6, 7, 8). 
 Each experimental unit consisted of five 1-week old seedlings grown together
 in a pot of nutrient solution during 14 days. 
 }


\usage{data(auxins)}

\format{
  A data frame with 150 observations on the following 5 variables.
  \describe{
    \item{\code{r}}{a numeric vector}
    \item{\code{h}}{a numeric vector}
    \item{\code{w}}{a numeric vector}
    \item{\code{y}}{a numeric vector}
    \item{\code{dose}}{a numeric vector}
  }
}

\details{
  The eight herbicide preparations are  naturally grouped into four pairs:
 (1, 5), (2, 6), (3, 7), and (4, 8), and in each pair the two herbicides
 should have the same active ingredients but different formulation
 constituents, which were assumed to be biologically inert. The data consist
 of the 150 observations y of dry weights, each observation being the weight
 of five plants grown in the same pot. All the eight herbicide preparations have
 essentially the same mode of action in the plant; they all act like the plant 
 auxins, which are plant regulators that affect cell enlongation an other 
 essential metabolic pathways. One of the objects of the experiment was to test
 if the response functions were identical except for a multiplicative factor in 
 the dose. This is a necessary, but not a sufficient, condition for a similar
  mode of action for the herbicides. The data are parts of a larger joint 
  action experiments with various herbicides and have also been used have been
   used in development of Random-Effect Models in Nonlinear Regression.

}

\source{
  Streibig, J. C. (1987). Joint action of root-absorbed mixtures of auxin 
  herbicides in Sinapis alba L. and barley (Hordeum vulgare L.). Weed Research 27, 337-347.

}

\references{
  Rudemo; M. Ruppert; D and Streibig, J. C. (1989). Random-Effect Models 
  in Nonlinear Regression with Applications to Bioassay.Biometrics, 45, 349-362.
}

\examples{
auxins.m1 <- drm(y ~ dose, h, adjust = "bc1", pmodels = data.frame(h,
+ h, 1, h), fct = LL.4(), data = auxins)
auxins.m2 <- drm(y ~ dose, h, adjust = "bc1", pmodels = data.frame(h,
+ 1, 1, h), fct = LL.4(), data = auxins)


}
\keyword{datasets}
