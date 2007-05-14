\name{drc-package}

\alias{drc-package}

\alias{drc}

\docType{package}

\title{

Modelling dose-response curves

}

\description{
Non-linear regression analysis for one or more curves with
focus on concentration-response, dose-response and
time-response curves used, for example, in biology, chemistry,
environmental sciences, medicine, pharmacology, toxicology.
}

\details{

\tabular{ll}{
Package: \tab drc\cr
LazyLoad: \tab yes\cr
LazyData: \tab yes\cr
Version: \tab 1.1-0\cr
License: \tab GNU GPL (version 2 or any later version)\cr
Depends: \tab R (>= 2.0.0), MASS\cr
Date: \tab 2007-02-01\cr
URL: \tab http://www.r-project.org, http://www.bioassay.dk\cr
Packaged: \tab Tue Aug 1 18:24:23 2006; ritz\cr
Built: \tab R 2.3.0; ; 2006-08-01 18:24:25; windows\cr
}

Index:
\preformatted{
ED                      Estimating effective dosage
FA                      Effect of Ferulic Acid on growth of ryegrass
MAX                     Maximum mean response
MM                      Michaelis-Menten model
PestSci                 Inhibition of photosynthesis
RScompetition           Competition between two biotypes
SI                      Comparing selectivity indices across curves
TM                      A suite of bioassays
TerMet                  The effect of terbuthylazin on growth rate
acidiq                  Acifluorfen and diquat tested on Lemna minor.
anova.drc               ANOVA for non-linear model fits
baro5                   The modified baro5 function
bcl3                    The three-parameter logistic function with
                        Brain-Cousens' modification
bcl4                    The four-parameter logistic function with
                        Brain-Cousens' modification
braincousens            The Brain-Cousens modification of the logistic
                        function
cabanne                 data name/kind
christensenetal2003     Potency of two herbicides
coef.drc                Extract Model Coefficients
colFct                  Collapse of levels
compParm                Comparison of parameters
confint.drc             Confidence Intervals for model parameters
daphnids                Daphnia test
diagnostics             Display information on estimation procedure
earthworms              Earthworms
finney71                Example from Finney (1971)
fitted.drc              Extract fitted values from model
glymet                  Glyphosate and metsulfuron-methyl tested on
                        algae.
heartrate               Heart rate baroreflexes for rabbits
hormesis                Hormesis data set
isobole                 Creating isobolograms
isomerization           Reaction rate of the catalytic isomerization
l2                      The two-parameter logistic function
l3                      The three-parameter logistic function
l4                      The four-parameter logistic function
l5                      The five-parameter logistic function
leaflength              Leaf length of barley
logLik.drc              Extracting the log likelihood
logistic                The logistic function
mdControl               Sets control arguments
mixdrc                  Fitting multiple dose-response curves with
                        random effects
mixture                 Fitting mixture models
ml3a                    The modified logistic function
ml3b                    The modified logistic function
ml3c                    The modified logistic function
ml4a                    The modified logistic function
ml4b                    The modified logistic function
ml4c                    The modified logistic function
mlogistic               The modified logistic function
multdrc                 Fitting multiple dose-response curves
onecurve                Root length
plot.drc                Plotting fitted curves for a drc object
plotraw                 Plotting raw dose response data.
predict.drc             Prediction
print.drc               Printing key features
print.summary.drc       Printing summary of non-linear model fits
residuals.drc           Extracting residuals from model
searchdrc               Searching through a range of initial parameter
                        values to obtain convergence
summary.drc             Summarising non-linear model fits
termec                  Mechlorprop and terbythylazine tested on Lemna
                        minor
twofactors              Dose-response curves for multiple factors
update.drc              Updating and re-fitting a model
vcov.drc                Calculating Variance-Covariance Matrix for
                        objects of class 'drc'
w2                      The two-parameter Weibull function
w3                      The three-parameter Weibull function
w4                      The four-parameter Weibull function
weibull                 The Weibull function
}

Further information is available in the following vignette:

\tabular{ll}{
\code{drc} \tab Bioassay analysis using R (source, pdf)\cr
}
}

\author{
  Christian Ritz <ritz@bioassay.dk>, 
  Jens Strebig <Jens.C.Streibig@agsci.kvl.dk>

  Maintainer: Christian Ritz <ritz@bioassay.dk>
}

\references{
  Ritz, C. and Streibig, J. C. (2005) Bioassay Analysis using R. \emph{J. Statist. Software}, \bold{12}, Issue 5.
}

\keyword{package}
\keyword{models}
\keyword{nonlinear}
%~~ Optionally other keywords from doc/KEYWORDS, one per line

\seealso{
%~~ Optional links to other man pages, e.g. ~~
%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
}

\examples{

model1 <- multdrc(rootl ~ conc, data = ryegrass)
plot(model1)

model2 <- multdrc(SLOPE ~ DOSE, CURVE, data = spinach)
plot(model2)

}
